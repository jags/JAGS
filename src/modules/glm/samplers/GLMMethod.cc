#include <config.h>
#include <rng/RNG.h>

#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include "GLMMethod.h"

#include <sampler/GraphView.h>
#include <sampler/Linear.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LinkNode.h>
#include <distribution/Distribution.h>
#include <rng/TruncatedNormal.h>
#include <module/ModuleError.h>

using std::string;
using std::vector;
using std::set;
using std::copy;
using std::sqrt;

extern cholmod_common *glm_wk;

static void getIndices(set<StochasticNode const *> const &schildren,
		       vector<StochasticNode const*> const &rows,
		       vector<int> &indices)
{
    indices.clear();
    
    for (unsigned int i = 0; i < rows.size(); ++i) {
	if (schildren.count(rows[i])) {
	    indices.push_back(i);
	}
    }
    
    if (indices.size() != schildren.size()) {
	throwLogicError("Size mismatch in getIndices");
    }
}

static Node const *getLinearPredictor(StochasticNode const *snode)
{
    Node const *lp = 0;

    switch(glm::GLMMethod::getFamily(snode)) {
    case GLM_NORMAL: case GLM_BERNOULLI: case GLM_BINOMIAL: case GLM_POISSON:
	lp = snode->parents()[0];
	break;
    case GLM_UNKNOWN:
	break;
    }
    
    LinkNode const *ln = dynamic_cast<LinkNode const*>(lp);
    if (ln) 
	lp = ln->parents()[0];
    
    return lp;
}

namespace glm {

    void GLMMethod::calDesign() const
    {
	vector<StochasticNode *> const &snodes = _view->nodes();
	vector<StochasticNode const *> const &schildren = 
	    _view->stochasticChildren();

	int *Xi = static_cast<int*>(_x->i);
	int *Xp = static_cast<int*>(_x->p);
	double *Xx = static_cast<double*>(_x->x);
	
	unsigned int nrow = schildren.size();
	unsigned int ncol = _view->length();
	if (nrow != _x->nrow || ncol != _x->ncol) {
	    throwLogicError("Dimension mismatch in GLMMethod::calDesign");
	}

	int c = 0; //column counter
	double *xnew = new double[_length_max];
	
	for (unsigned int i = 0; i < snodes.size(); ++i) {

	    unsigned int length = snodes[i]->length();

	    if (_init || !_fixed[i]) {

		for (unsigned int j = 0; j < length; ++j) {
		    for (int r = Xp[c+j]; r < Xp[c+j+1]; ++r) {
			Xx[r] = -getMean(Xi[r]);
		    }
		}
		
		double const *xold = snodes[i]->value(_chain);	    
		copy(xold, xold + length, xnew);
		
		for (unsigned int j = 0; j < length; ++j) {
		    
		    xnew[j] += 1;
		    _sub_views[i]->setValue(xnew, length, _chain);
		    for (int r = Xp[c+j]; r < Xp[c+j+1]; ++r) {
			Xx[r] += getMean(Xi[r]);
		    }
		    xnew[j] -= 1;
		}
		_sub_views[i]->setValue(xnew, length, _chain);
	    }
	    
	    c += length;
	}

	delete [] xnew;
	
    }
    
    GLMMethod::GLMMethod(GraphView const *view, 
			 vector<GraphView const *> const &sub_views,
			 unsigned int chain, bool link)
	: _lp(view->stochasticChildren().size()),
	  _view(view), _chain(chain), _sub_views(sub_views),
	  _x(0), _factor(0), _fixed(sub_views.size(), false), 
	  _length_max(0), _nz_prior(0), _init(true)
    {
	vector<StochasticNode const*> const &schildren = 
	    view->stochasticChildren();

	int nrow = schildren.size();
	int ncol = view->length();

	//Set up linear predictor
	for (int i = 0; i < nrow; ++i) {
	    _lp[i] = getLinearPredictor(schildren[i])->value(chain);
	}

	vector<int> Xp(ncol + 1);
	vector<int> Xi;
    
	int c = 0; //column counter
	int r = 0; //count of number of non-zero entries

	for (unsigned int p = 0; p < _sub_views.size(); ++p) {

	    set<StochasticNode const *> children_p;
	    children_p.insert(sub_views[p]->stochasticChildren().begin(),
			      sub_views[p]->stochasticChildren().end());
	    vector<int> indices;
	    getIndices(children_p, schildren, indices);

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
	_x = cholmod_allocate_sparse(nrow, ncol, r, 1, 1, 0, CHOLMOD_REAL, glm_wk);
	int *_xp = static_cast<int*>(_x->p);
	int *_xi = static_cast<int*>(_x->i);

	copy(Xp.begin(), Xp.end(), _xp);
	copy(Xi.begin(), Xi.end(), _xi);

	// Check for constant linear terms
	for (unsigned int i = 0; i < sub_views.size(); ++i) {
	    _fixed[i] = checkLinear(sub_views[i], true, link);
	}
    }

    GLMMethod::~GLMMethod()
    {
	cholmod_free_sparse(&_x, glm_wk);
    }
    
    /* 
       Symbolic analysis of the posterior precision matrix for the
       Cholesky decomposition.
       
       This only needs to be done once, when the GLMMethod is
       craeted. It is a stripped-down version of the code in update.
       Note that the values of the sparse matrices are never
       referenced.
    */
    void GLMMethod::symbolic()  
    {
	unsigned int nrow = _view->length();

	// Prior contribution 
	cholmod_sparse *Aprior =  
	    cholmod_allocate_sparse(nrow, nrow, _nz_prior, 1, 1, 0, CHOLMOD_PATTERN, glm_wk); 
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

    void GLMMethod::calCoef(double *&b, cholmod_sparse *&A, double Temp)
    {
	//   The log of the full conditional density takes the form
	//   -(t(x) %*% A %*% x - 2 * b %*% x)/2
	//   where A is the posterior precision and the mean mu solves
	//   A %*% mu = b

	//   For computational convenience we take xold, the current value
	//   of the sampled nodes, as the origin

	unsigned int nrow = _view->length();
	b = new double[nrow];

	cholmod_sparse *Aprior =  
	    cholmod_allocate_sparse(nrow, nrow, _nz_prior, 1, 1, 0, 
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

	cholmod_sparse *t_x = cholmod_transpose(_x, 1, glm_wk);
	int *Tp = static_cast<int*>(t_x->p);
	int *Ti = static_cast<int*>(t_x->i);
	double *Tx = static_cast<double*>(t_x->x);

	for (unsigned int c = 0; c < t_x->ncol; ++c) {
	    double tau = Temp * getPrecision(c);
	    double delta = tau * (getValue(c) - getMean(c));
	    double sigma = sqrt(tau);
	    for (int r = Tp[c]; r < Tp[c+1]; ++r) {
		b[Ti[r]] += Tx[r] * delta;
		Tx[r] *= sigma;
	    }
	}

	cholmod_sparse *Alik = cholmod_aat(t_x, 0, 0, 1, glm_wk);
	cholmod_free_sparse(&t_x, glm_wk);
	double alpha[2] = {1, 0};
	double beta[2] = {1, 0};
	
	A = cholmod_add(Aprior, Alik, one, one, 1, 0, glm_wk);
	cholmod_free_sparse(&Aprior, glm_wk);
	cholmod_free_sparse(&Alik, glm_wk);
    }

    void GLMMethod::updateLM(RNG *rng, double Temp) 
    {
	//   The log of the full conditional density takes the form
	//   -(t(x) %*% A %*% x - 2 * b %*% x)/2
	//   where A is the posterior precision and the mean mu solves
	//   A %*% mu = b

	//   For computational convenience we take xold, the current value
	//   of the sampled nodes, as the origin

	if (_init) {
	    calDesign();
	    symbolic();
	    _init = false;
	}

	double *b = 0;
	cholmod_sparse *A = 0;
	calCoef(b, A, Temp);
	
	// Get LDL' decomposition of posterior precision
	A->stype = -1;
	int ok = cholmod_factorize(A, _factor, glm_wk);
	cholmod_free_sparse(&A, glm_wk);
	if (!ok) {
	    throwRuntimeError("Cholesky decomposition failure in GLMMethod");
	}

	// Use the LDL' decomposition to generate a new sample
	// with mean mu such that A %*% mu = b and precision A. 
	
	unsigned int nrow = _view->length();
	cholmod_dense *w = cholmod_allocate_dense(nrow, 1, nrow, CHOLMOD_REAL, 
						  glm_wk);

	// Permute RHS
	double *wx = static_cast<double*>(w->x);
	int *perm = static_cast<int*>(_factor->Perm);
	for (unsigned int i = 0; i < nrow; ++i) {
	    wx[i] = b[perm[i]];
	}

	cholmod_dense *u1 = cholmod_solve(CHOLMOD_L, _factor, w, glm_wk);
	updateAuxiliary(u1, _factor, rng);

	double *u1x = static_cast<double*>(u1->x);
	if (_factor->is_ll) {
	    // LL' decomposition
	    for (unsigned int r = 0; r < nrow; ++r) {
		u1x[r] += rng->normal();
	    }
	}
	else {
	    // LDL' decomposition. The diagonal D matrix is stored
	    // as the diagonal of _factor
	    int *fp = static_cast<int*>(_factor->p);
	    double *fx = static_cast<double*>(_factor->x);
	    for (unsigned int r = 0; r < nrow; ++r) {
		u1x[r] += rng->normal() * sqrt(fx[fp[r]]);
	    }
	}

	cholmod_dense *u2 = cholmod_solve(CHOLMOD_DLt, _factor, u1, glm_wk);

	// Permute solution
	double *u2x = static_cast<double*>(u2->x);
	for (unsigned int i = 0; i < nrow; ++i) {
	    b[perm[i]] = u2x[i];
	}

	cholmod_free_dense(&w, glm_wk);
	cholmod_free_dense(&u1, glm_wk);
	cholmod_free_dense(&u2, glm_wk);

	//Shift origin back to original scale
	int r = 0;
	for (vector<StochasticNode*>::const_iterator p = 
		 _view->nodes().begin();  p != _view->nodes().end(); ++p)
	{
	    unsigned int length = (*p)->length();
	    double const *xold = (*p)->value(_chain);
	    for (unsigned int i = 0; i < length; ++i, ++r) {
		b[r] += xold[i];
	    }
	}

	_view->setValue(b, nrow, _chain);
	delete [] b;
    }

    void GLMMethod::updateLMGibbs(RNG *rng) 
    {
	// Update element-wise. Less efficient than updateLM but
	// does not require a Cholesky decomposition, and is 
	// necessary for truncated parameters

	if (_init) {
	    if (_view->length() != _sub_views.size()) {
		throwLogicError("updateLMGibbs can only act on scalar nodes");
	    }
	    calDesign();
	    _init = false;
	}

	double *b = 0;
	cholmod_sparse *A = 0;
	calCoef(b, A);

	int nrow = _view->length();
	vector<double> theta(nrow);
	_view->getValue(theta, _chain);

	int *Ap = static_cast<int*>(A->p);
	int *Ai = static_cast<int*>(A->i);
	double *Ax = static_cast<double*>(A->x);

	//Extract diagonal from A
	vector<double> diagA(nrow);
	for (int c = 0; c < nrow; ++c) {
	    for (int j = Ap[c]; j < Ap[c+1]; ++j) {
		if (Ai[j] == c) {
		    diagA[c] = Ax[j];
		    break;
		}
	    }
	}

	//Update element-wise
	for (int i = 0; i < nrow; ++i) {
	    
	    double theta_old = theta[i];
		
	    double mu  = theta[i] + b[i]/diagA[i];
	    double sigma = sqrt(1/diagA[i]);
	    StochasticNode const *snode = _sub_views[i]->nodes()[0];
	    double const *l = snode->lowerLimit(_chain);
	    double const *u = snode->upperLimit(_chain);
		
	    if (l && u) {
		theta[i] = inormal(*l, *u, rng, mu, sigma);
	    }
	    else if (l) {
		theta[i] = lnormal(*l, rng, mu, sigma);
	    }
	    else if (u) {
		theta[i] = rnormal(*u, rng, mu, sigma);
	    }
	    else {
		theta[i] = mu + rng->normal() * sigma;
	    }
		
	    double delta = theta[i] - theta_old;
	    for (int j = Ap[i]; j < Ap[i+1]; ++j) {
		b[Ai[j]] -= delta * Ax[j];
	    }
	}

	_view->setValue(theta,  _chain);
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

    double GLMMethod::getMean(unsigned int i) const
    {
	return *_lp[i];
    }

    GLMFamily GLMMethod::getFamily(StochasticNode const *snode)
    {
	string const &name = snode->distribution()->name();
	if (name == "dbern") {
	    return GLM_BERNOULLI;
	}
	else if (name == "dbin") {
	    return GLM_BINOMIAL;
	}
	else if (name == "dpois") {
	    return GLM_POISSON;
	}
	else if (name == "dnorm") {
	    return GLM_NORMAL;
	}
	else {
	    return GLM_UNKNOWN;
	}
    }

    void GLMMethod::updateAuxiliary(cholmod_dense *b, cholmod_factor *N, RNG *rng)
    {
    }

}
