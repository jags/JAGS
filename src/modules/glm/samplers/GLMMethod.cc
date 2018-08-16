#include <config.h>

#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include "GLMMethod.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <sampler/Linear.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
//#include <graph/LinkNode.h>
#include <distribution/Distribution.h>
#include <rng/TruncatedNormal.h>
#include <module/ModuleError.h>
#include <util/logical.h>
#include <rng/RNG.h>

using std::string;
using std::vector;
using std::set;
using std::copy;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {

static void getIndices(set<StochasticNode *> const &schildren,
		       vector<StochasticNode *> const &allchildren,
		       vector<unsigned int > const &allrows,
		       vector<int> &indices)
{
    indices.clear();
    
    for (unsigned int i = 0; i < allchildren.size(); ++i) {
	if (schildren.count(allchildren[i])) {
	    for (unsigned int j = allrows[i]; j < allrows[i+1]; ++j) {
		indices.push_back(j);
	    }
	}
    }
}

namespace glm {

    void GLMMethod::calDesign() const
    {
	if (allTrue(_fixed)) return; //Move along, nothing to see here
	
	vector<StochasticNode *> const &snodes = _view->nodes();

	int *Xi = static_cast<int*>(_x->i);
	int *Xp = static_cast<int*>(_x->p);
	double *Xx = static_cast<double*>(_x->x);
	
	unsigned int nrow = _x->nrow;
	unsigned int ncol = _view->length();
	if (ncol != _x->ncol) {
	    throwLogicError("Dimension mismatch in GLMMethod::calDesign");
	}

	int c = 0; //column counter
	vector<double> xnew(_length_max);

	//Need to set up these vectors so we can map back from a row of
	//the design matrix to the corresponding outcome. Note that this
	//is only need when we have multivariate outcomes; with only
	//scalar outcomes the mapping is trivial
	vector<double const*> outcome_ptr(nrow);
	vector<unsigned int> outcome_idx(nrow);
	unsigned int r = 0; //row counter
	for (unsigned int i = 0; i < _outcomes.size(); ++i) {
	    double const *vm = _outcomes[i]->vmean();
	    for (unsigned int j = 0; j < _outcomes[i]->length(); ++j, ++r) {
		outcome_ptr[r] = vm;
		outcome_idx[r] = j;
	    }
	}
 	
	for (unsigned int i = 0; i < snodes.size(); ++i) {

	    unsigned int length = snodes[i]->length();

	    if (!_fixed[i]) {

		for (unsigned int j = 0; j < length; ++j) {
		    for (int r = Xp[c+j]; r < Xp[c+j+1]; ++r) {
 			unsigned int row = Xi[r];
			Xx[r] = - outcome_ptr[row][outcome_idx[row]];
		    }
		}
		
		double const *xold = snodes[i]->value(_chain);	    
		copy(xold, xold + length, xnew.begin());
		
		for (unsigned int j = 0; j < length; ++j) {
		    double xsave = xnew[j];
		    xnew[j] += 1;
		    _sub_views[i]->setValue(&xnew[0], length, _chain);
		    for (int r = Xp[c+j]; r < Xp[c+j+1]; ++r) {
			unsigned int row = Xi[r];
			Xx[r] += outcome_ptr[row][outcome_idx[row]];
		    }
		    xnew[j] = xsave;
		}
		_sub_views[i]->setValue(&xnew[0], length, _chain);
	    }
	    
	    c += length;
	}
    }
    
    GLMMethod::GLMMethod(GraphView const *view, 
			 vector<SingletonGraphView const *> const &sub_views,
			 vector<Outcome *> const &outcomes,
			 unsigned int chain)
	: _view(view), _chain(chain), _sub_views(sub_views),
	  _outcomes(outcomes),
	  _x(0), _factor(0), _fixed(sub_views.size(), false), 
	  _length_max(0), _nz_prior(0)
    {
	view->checkFinite(chain); //Check validity of initial values
	
	vector<StochasticNode *> const &schildren = 
	    view->stochasticChildren();

	vector<unsigned int> rows(schildren.size() + 1);
	rows[0] = 0;
	for (unsigned int i = 0; i < schildren.size(); ++i) {
	    rows[i+1] = rows[i] + schildren[i]->length();
	}
	int ncol = view->length();
	int nrow = rows[schildren.size()];
	
	vector<int> Xp(ncol + 1);
	vector<int> Xi;
    
	int c = 0; //column counter
	int r = 0; //count of number of non-zero entries

	for (unsigned int p = 0; p < _sub_views.size(); ++p) {

	    set<StochasticNode *> children_p;
	    children_p.insert(sub_views[p]->stochasticChildren().begin(),
			      sub_views[p]->stochasticChildren().end());
	    vector<int> indices;
	    getIndices(children_p, schildren, rows, indices);

	    unsigned int length = _sub_views[p]->length();
	    for (unsigned int i = 0; i < length; ++i, ++c) {
		Xp[c] = r;
		for (unsigned int j = 0; j < indices.size(); ++j, ++r) {
		    Xi.push_back(indices[j]);
		}
	    }

	    //Save these values for later calculations
	    _nz_prior += length * length; //No. of non-zeros in prior precision
	    if (length > _length_max) {
		_length_max = length; //Length of longest sampled node
	    }
	}
	Xp[c] = r;

	//Set up sparse representation of the design matrix
	#pragma omp critical
	_x = cholmod_allocate_sparse(nrow, ncol, r, 1, 1, 0, CHOLMOD_REAL, glm_wk);
	int *_xp = static_cast<int*>(_x->p);
	int *_xi = static_cast<int*>(_x->i);

	copy(Xp.begin(), Xp.end(), _xp);
	copy(Xi.begin(), Xi.end(), _xi);

	// At this point, all elements of _fixed are set to false, so
	// a call to calDesign calculates the whole design matrix
	calDesign();
	
	// In future calls to calDesign, we do not want to recalculate
	// fixed linear terms.
	for (unsigned int i = 0; i < sub_views.size(); ++i) {
	    // FIXME: For future reference, we will need to make sure this
	    // still works correctly for log-linear models.
	    _fixed[i] = checkLinear(sub_views[i], true, true);
	}
    }

    GLMMethod::~GLMMethod()
    {
	while(!_outcomes.empty()) {
	    delete _outcomes.back();
	    _outcomes.pop_back();
	}
	#pragma omp critical
	cholmod_free_sparse(&_x, glm_wk);
    }
    
    /* 
       Symbolic analysis of the posterior precision matrix for the
       Cholesky decomposition.
       
       This only needs to be done once, when the GLMMethod is
       created. It is a stripped-down version of the code in update.
       Note that the values of the sparse matrices are never
       referenced.
    */
    void GLMMethod::symbolic()  
    {
	unsigned int nrow = _view->length();

	// Prior contribution 
	cholmod_sparse *Aprior = cholmod_allocate_sparse(nrow, nrow, _nz_prior, 1, 1, 0, CHOLMOD_PATTERN, glm_wk); 
	int *Ap = static_cast<int*>(Aprior->p);
	int *Ai = static_cast<int*>(Aprior->i);

	int c = 0;
	int r = 0;
	vector<StochasticNode*> const &snodes = _view->nodes();
	for (vector<StochasticNode*>::const_iterator p = snodes.begin();
	     p != snodes.end(); ++p)
	{
	    StochasticNode *snode = *p;
	    unsigned int length = snode->length();

	    /* 
	       Fixme: we're assuming the prior precision of each node
	       is dense, whereas it may be sparse.
	    */
	    int cbase = c; //first column in this diagonal block
	    for (unsigned int j = 0; j < length; ++j, ++c) {
		Ap[c] = r;
		for (unsigned int i = 0; i < length; ++i, ++r) {
		    Ai[r] = cbase + i;
		}
	    }
	}
	Ap[c] = r;
	
	// Likelihood contribution

	cholmod_sparse *t_x = cholmod_transpose(_x, 0, glm_wk);
	cholmod_sort(t_x, glm_wk);
	cholmod_sparse *Alik = cholmod_aat(t_x, 0, 0, 0, glm_wk);
	cholmod_sparse *A = cholmod_add(Aprior, Alik, 0, 0, 0, 0, glm_wk);
	
	//Free working matrices
	cholmod_free_sparse(&t_x, glm_wk);
	cholmod_free_sparse(&Aprior, glm_wk);
	cholmod_free_sparse(&Alik, glm_wk);
	
	A->stype = -1;
	_factor = cholmod_analyze(A, glm_wk); 
	cholmod_free_sparse(&A, glm_wk);
    }

    void GLMMethod::calCoef(double *&b, cholmod_sparse *&A) 
    {
	//   The log of the full conditional density takes the form
	//   -(t(x) %*% A %*% x - 2 * b %*% x)/2
	//   where A is the posterior precision and the mean mu solves
	//   A %*% mu = b

	//   For computational convenience we take xold, the current value
	//   of the sampled nodes, as the origin

	unsigned int nrow = _view->length();
	b = new double[nrow];

	cholmod_sparse *Aprior = 0;
	#pragma omp critical
	Aprior = cholmod_allocate_sparse(nrow, nrow, _nz_prior, 1, 1, 0,
					 CHOLMOD_REAL, glm_wk); 
    
	// Set up prior contributions to A, b
	int *Ap = static_cast<int*>(Aprior->p);
	int *Ai = static_cast<int*>(Aprior->i);
	double *Ax = static_cast<double*>(Aprior->x);

	//FIXME. We are assuming contributions to prior are dense
	int c = 0;
	int r = 0;
	vector<StochasticNode*> const &snodes = _view->nodes();
	for (vector<StochasticNode*>::const_iterator p = snodes.begin();
	     p != snodes.end(); ++p)
	{
	    StochasticNode *snode = *p;
	    double const *priormean = snode->parents()[0]->value(_chain);
	    double const *priorprec = snode->parents()[1]->value(_chain);
	    double const *xold = snode->value(_chain);
	    unsigned int length = snode->length();
	
	    int cbase = c; //first column of this diagonal block
	    for (unsigned int i = 0; i < length; ++i, ++c) {
		b[c] = 0;
		Ap[c] = r;
		for (unsigned int j = 0; j < length; ++j, ++r) {
		    b[c] += priorprec[i + length*j] * (priormean[j] - xold[j]);
		    Ai[r] = cbase + j;
		    Ax[r] = priorprec[i + length*j];
		}
	    }
	}
	Ap[c] = r;

	// Recalculate the design matrix, if necessary
	calDesign();

	// Likelihood contributions
	//   
	//   b += t(X) %*% tau %*% (Y - mu)
	//   A += t(X) %*% tau %*% X
    
	//   where 
	//   - X is the design matrix
	//   - tau is the (diagonal) variance matrix of the stochastic children
	//   - mu is the mean of the stochastic children
	//   - Y is the value of the stochastic children

	/* FIXME: We can also transpose _x into a pre-defined sparse
	   matrix. Try this with the lsat example to see what speed-up
	   it gives.
	*/
	cholmod_sparse *t_x = 0;
	#pragma omp critical
	{
	    t_x = cholmod_transpose(_x, 1, glm_wk);
	    cholmod_sort(t_x, glm_wk); //Needed for multivariate outcomes
	}
	
	int *Tp = static_cast<int*>(t_x->p);
	int *Ti = static_cast<int*>(t_x->i);
	double *Tx = static_cast<double*>(t_x->x);

	/*
	//FIXME: reusing previously defined c,r here
	for (c = 0; c < t_x->ncol; ++c) {
	    double tau = _outcomes[c]->precision();
	    double delta = tau * (_outcomes[c]->value() - _outcomes[c]->mean());
	    //double sigma = sqrt(tau);
	    for (r = Tp[c]; r < Tp[c+1]; ++r) {
		b[Ti[r]] += Tx[r] * delta;
		//Tx[r] *= sigma;
		Tx[r] *= tau;
	    }
	}
	*/
	
	//FIXME: reusing previously defined c,r here
	c = 0;
	for (unsigned int i = 0; i < _outcomes.size(); ++i) {
	    unsigned int m = _outcomes[i]->length();
	    if (m == 1) {
		//Scalar outcome
		double tau = _outcomes[i]->precision();
		double delta = tau * (_outcomes[i]->value() -
				      _outcomes[i]->mean());
		for (r = Tp[c]; r < Tp[c+1]; ++r) {
		    b[Ti[r]] += Tx[r] * delta;
		    Tx[r] *= tau;
		}
	    }
	    else {
		//Multivariate outcome
		double const *tau = _outcomes[i]->vprecision();
		double const *Y = _outcomes[i]->vvalue();
		double const *mu = _outcomes[i]->vmean();

		vector<double> delta(m,0);
		for (unsigned int j = 0; j < m; ++j) {
		    for (unsigned int k = 0; k < m; ++k) {
			delta[j] += tau[m*j+k] * (Y[k] - mu[k]);
		    }
		}

		/* Sanity checks */
		int nzrow = Tp[c+1] - Tp[c];
		for (r = Tp[c]; r < Tp[c+1]; ++r) {
		    int row = Ti[r];
		    for (unsigned int j = 0; j < m; ++j) {
			unsigned int col = c + j;
			if (Tp[col+1] - Tp[col] != nzrow) {
			    string ms("non-zero pattern mismatch in GLMMethod");
			    throwLogicError(ms);
			}
			if (Ti[r+j*nzrow] != row) {
			    string ms("row order mismatch in GLMMethod");
			    throwLogicError(ms);
			}
		    }
		}

		/* We assume that all columns of t_x corresponding to
		   the same outcome have the same non-zero pattern
		   and the columns are sorted.
		*/
		vector<double> TxTau(m); // t_x %*% tau
		for (r = Tp[c]; r < Tp[c+1]; ++r) {
		    for (unsigned int j = 0; j < m; ++j) {
			b[Ti[r]] += Tx[r + j*nzrow] * delta[j];
			TxTau[j] = 0.0;
			for (unsigned int k = 0; k < m; ++k) {
			    TxTau[j] += Tx[r + k*nzrow] * tau[m*j+k];
			}
		    }
		    
		    for (unsigned int j = 0; j < m; ++j) {
			Tx[r + j*nzrow] = TxTau[j];
		    }
		}
	    }
	    c += m;
	}

	#pragma omp critical
	{
	    cholmod_sparse *Alik = cholmod_ssmult(t_x, _x, CHOLMOD_REAL, 1, 0,
						  glm_wk);
	    cholmod_free_sparse(&t_x, glm_wk);
	    double one[2] = {1, 0};
	    A = cholmod_add(Aprior, Alik, one, one, 1, 0, glm_wk);
	    
	    cholmod_free_sparse(&Aprior, glm_wk);
	    cholmod_free_sparse(&Alik, glm_wk);
	}
    }

    bool GLMMethod::isAdaptive() const
    {
	return false;
    }
    
    void GLMMethod::adaptOff()
    {
    }

    bool GLMMethod::checkAdaptation() const
    {
	return true;
    }

}}
