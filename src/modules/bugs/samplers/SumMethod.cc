#include <config.h>
#include <rng/RNG.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/Graph.h>
#include <sampler/Linear.h>
#include <util/nainf.h>
#include <graph/NodeError.h>

#include "SumMethod.h"
#include <cmath>
#include <numeric>
#include <cfloat>
#include <stdexcept>
#include <set>
#include <algorithm>

using std::vector;
using std::accumulate;
using std::min;
using std::max;
using std::ceil;
using std::logic_error;
using std::floor;
using std::fabs;
using std::list;
using std::find;

namespace jags {
    
    namespace bugs {

	template<class T>
	double sumValue(vector<T*> const &x, unsigned int ch)
	{
	    double y = 0;
	    for (typename vector<T*>::const_iterator p = x.begin();
		 p != x.end(); ++p)
	    {
		double const *v = (*p)->value(ch);
		unsigned long n = (*p)->length();
		y = accumulate(v, v+n, y);
	    }
	    return y;
	}

	/*
	//Test to see of two sets intersect
	template <typename T, typename U>
	bool intersects(set<T*> const &st, set<U*> const &su)
        {
	    typename set<T*>::const_iterator t = st.begin();
	    typename set<U*>::const_iterator u = su.begin();

	    while (t != st.end() && u != su.end()) {
		if (*t < *u) {
		    if (t != st.end()) ++t;
		}
		else if (*t > *u) {
		    if (u != su.end()) ++u;
		}
		else {
		    return true;
		}
	    }
	    return false;
	}
	*/
	    
	StochasticNode *
	SumMethod::isCandidate(StochasticNode *snode, Graph const &graph)
	{
	    //We can only sample scalar nodes
	    if (snode->length() != 1) return 0;
	    
	    SingletonGraphView gv(snode, graph);

	    /* 
	       Check stochastic children
	       
	    */
	    vector<StochasticNode*> const &schildren = gv.stochasticChildren();
	    StochasticNode *sumchild = 0;
	    for (unsigned int i = 0; i < schildren.size(); ++i) {
		if (schildren[i]->distribution()->name() == "sum") {
		    /* There should be only one sum node and it must be
		       observed */
		    if (sumchild) return 0;
		    if(!isObserved(schildren[i])) return 0;
		    sumchild = schildren[i];
		}
	    }
	    if (!sumchild) return 0;
	    
	    /* Deterministic descendants must be an additive function
	     * of snode
	     */
	    if (schildren.size() == 1) {
		if (!checkAdditive(&gv, false)) return 0;
	    }
	    else {
		/* 
		   We only care that sumnode is an additive function
		   so construct a local graph that excludes other
		   stochastic children for additivity check
		*/
		Graph lgraph;
		lgraph.insert(snode);
		lgraph.insert(sumchild);
		vector<DeterministicNode*> const &dchildren =
		    gv.deterministicChildren();
		for (unsigned int j = 0; j < dchildren.size(); ++j) {
		    lgraph.insert(dchildren[j]);
		}
		SingletonGraphView lgv(snode, lgraph);
		if (!checkAdditive(&lgv, false)) return 0;
	    }

	    return sumchild;
	}
	
	bool SumMethod::canSample(vector<StochasticNode*> const &snodes,
				  Graph const &graph)
	{
	    //Are individual nodes candidates?
	    StochasticNode *sumchild = isCandidate(snodes[0], graph);
	    if (sumchild == 0) return false;
	    for (unsigned int i = 1; i < snodes.size(); ++i) {
		if (isCandidate(snodes[i], graph) != sumchild) return false;
	    }

	    //Check consistency of discreteness
	    bool discrete = sumchild->isDiscreteValued();
	    for (unsigned int i = 1; i < snodes.size(); ++i) {
		if (snodes[i]->isDiscreteValued() != discrete) return false;
	    }

	    //Together are the nodes additive with a fixed intercept?

	    GraphView gv(snodes, graph, true);
	    
	    /* Only deterministic nodes with a path to sumnode need to
	       be additive. Constract a local graph that strips out
	       the rest.
	    */
	    Graph lgraph; //local graph
	    vector<DeterministicNode*> const &dc = gv.deterministicChildren();
	    for (vector<DeterministicNode*>::const_reverse_iterator p =
		     dc.rbegin(); p != dc.rend(); ++p)
	    {
		list<StochasticNode*> const *psc = (*p)->stochasticChildren();
		if (find(psc->begin(), psc->end(), sumchild) != psc->end()) {
		    lgraph.insert(*p);
		}
		else {
		    list<DeterministicNode*> const *pdc =
			(*p)->deterministicChildren();
		    for(list<DeterministicNode*>::const_iterator q =
			    pdc->begin(); q != pdc->end(); ++q)
		    {
			if (lgraph.contains(*q)) {
			    lgraph.insert(*p);
			    break;
			}
		    }
		}
	    }
	    lgraph.insert(sumchild);
	    for (unsigned int i = 0; i < snodes.size(); ++i) {
		lgraph.insert(snodes[i]);
	    }

	    return checkAdditive(snodes, lgraph, true);
	}

	SumMethod::SumMethod(GraphView const *gv, unsigned int chain)
	    : MutableSampleMethod(), _gv(gv), _chain(chain),
	      _sum(gv->stochasticChildren()[0]->value(chain)[0]),
	      _discrete(gv->nodes()[0]->isDiscreteValued()),
	      _x(gv->length()), _i(0), _j(0), _sumchild(0), _fast(false),
	      _sumdiff(0), _iter(0), _width(2), _max(10), _adapt(true)
	{
	    if (gv->stochasticChildren().size() == 1) {
		_sumchild = gv->stochasticChildren().front();
		_fast = true;
	    }
	    else {
		vector<StochasticNode*> const &schildren =
		    gv->stochasticChildren();
		for (unsigned int i = 0; i < schildren.size(); ++i) {
		    if (schildren[i]->distribution()->name() == "sum") {
			_sumchild = schildren[i];
			break;
		    }
		}
	    }
	    
	    gv->getValue(_x, chain);
	    
	    if (_sumchild->logDensity(chain, PDF_LIKELIHOOD) != 0) {
		//If initial values are inconsistent with outcome we
		//try to adjust them

		// Calculate intercept (usually zero)
		double y = sumValue<const Node>(_sumchild->parents(), chain);
		double sumx = sumValue<StochasticNode>(gv->nodes(), chain);
		double alpha = y - sumx;

		// Calculate target sum for sampled values
		double sumx_new = _sumchild->value(chain)[0] - alpha;

		// Rescale
		unsigned long N = _x.size();
		vector<double> xnew;
		if (_discrete) {
		    xnew = vector<double>(N, floor(sumx_new/N));
		    double delta = sumx_new -
			accumulate(xnew.begin(), xnew.end(), 0);
		    xnew[N-1] += delta;
		}
		else {
		    xnew = vector<double>(N, sumx_new/N);
		}

		gv->setValue(xnew, chain);
		if (_sumchild->logDensity(chain, PDF_LIKELIHOOD) != 0) {
		    throw logic_error("SumMethod failed to fix initial values");
		}
		if (jags_finite(gv->logFullConditional(_chain))) {
		    _x = xnew; //Preserve changes
		}
		else {
		    throw NodeError(_sumchild, "SumMethod cannot fix the stochastic parents of this node to satisfy the sum constraint.\nYou must supply initial values for the parents");
		}
	    }

	    //Check validity of initial values
	    gv->checkFinite(chain);
	}
	
	SumMethod::~SumMethod()
	{
	}

	void SumMethod::updateStep(RNG *rng)
	{
	    // Test current value
	    double g0 = logDensity();
	    if (!jags_finite(g0)) {
		return;
	    }

	    // Generate auxiliary variable
	    double z = g0 - rng->exponential();

	    //Discrete width
	    double width = _discrete ? ceil(_width) : _width;
	    
	    // Generate random interval of width "_width" about current value
	    double left = _discrete ? floor(rng->uniform() * (width + 1)) :
		rng->uniform() * width;
	    double xold = value();
	    double L = xold - left;
	    double R = L + width;

	    double lower = JAGS_NEGINF, upper = JAGS_POSINF;
	    getLimits(&lower, &upper);
	    if (lower == upper) {
		// We can't move this pair, so skip them
		return;
	    }
		    
	    // Stepping out 

	    // Randomly set number of steps in left and right directions,
	    // subject to the limit in the maximal size of the interval
	    int j = static_cast<int>(rng->uniform() * _max);
	    int k = _max - 1 - j;

	    if (L < lower) {
		L = lower;
	    }
	    else {
		setValue(L);
		while (j-- > 0 && logDensity() > z) {
		    L -= width;
		    if (L < lower) {
			L = lower;
			break;
		    }
		    setValue(L);
		}
	    }

	    if (R > upper) {
		R = upper;
	    }
	    else {
		setValue(R);
		while (k-- > 0 && logDensity() > z) {
		    R += width;
		    if (R > upper) {
			R = upper;
			break;
		    }
		    setValue(R);
		}
	    }

	    // Keep sampling from the interval until acceptance (the loop is
	    // guaranteed to terminate).
	    double xnew;
	    for(;;) {
		if (_discrete) {
		    xnew =  L + floor(rng->uniform() * (R - L + 1));
		}
		else {
		    xnew =  L + rng->uniform() * (R - L);
		}
		setValue(xnew);
		double g = logDensity();
		if (g >= z - DBL_EPSILON) {
		    // Accept point
		    break;
		}
		else {
		    // shrink the interval
		    if (xnew < xold) {
			L = xnew;
		    }
		    else {
			R = xnew;
		    }
		}
	    }
	    
	    if (_adapt) {
		_sumdiff += fabs(xnew - xold);
	    }
	}

	void SumMethod::update(RNG *rng)
	{
	    unsigned int len = _gv->length();

	    if (len == 1) {
		//Value of single sampled node is determined by sum
		//constraint so it can never move. But we still check
		//that the sum constraint is satisfied.
		if (_sumchild->logDensity(_chain, PDF_LIKELIHOOD) != 0) {
		    throw logic_error("Failure to preserve sum in SumMethod");
		}
		return;
	    }

	    //Generate random permutation (rp) of indices 0...(len-1)
	    //using Fisher Yates algorithm.
	    vector<unsigned int> rp(len+1);
	    for (unsigned int i = 0; i < len; ++i) {
		unsigned int j =
		    static_cast<unsigned int>(rng->uniform() * (i+1));
		rp[i] = rp[j];
		rp[j] = i;
	    }
	    rp[len] = rp[0]; //Cycle back to the start

	    for (unsigned int i = 0; i < len; ++i) {
		_i = rp[i];
		_j = rp[i+1];
		updateStep(rng); //Slice update of par (_i, _j)
	    }
	    
	    if (_adapt) {
		if (++_iter % 50 == 0) {
		    _width = _sumdiff / (50 * len);
		    _sumdiff = 0;
		    if (_discrete) {
			_width = ceil(_width);
		    }
		}
	    }

	    //updateStep does not update the deterministic descendants
	    //but we need to leave these in a correct state in case they
	    //are being monitored
	    vector<DeterministicNode*> const &d = _gv->deterministicChildren();
	    for(vector<DeterministicNode*>::const_iterator p = d.begin();
		p != d.end(); ++p)
	    {
		(*p)->deterministicSample(_chain);
	    }

	    //Sanity check
	    if (_sumchild->logDensity(_chain, PDF_LIKELIHOOD) != 0) {
		throw logic_error("Failure to preserve sum in SumMethod");
	    }

	}

	void SumMethod::setValue(double x)
	{
	    double delta = x - _x[_i];
	    _x[_i] = x;      
	    _x[_j] -= delta; 

	    _gv->nodes()[_i]->setValue(&_x[_i], 1, _chain);
	    _gv->nodes()[_j]->setValue(&_x[_j], 1, _chain);

	    if (!_fast) {
		//Have to update deterministic descendants
		vector<DeterministicNode*> const &dtrm =
		    _gv->deterministicChildren();
		for (vector<DeterministicNode*>::const_iterator p(dtrm.begin());
		     p != dtrm.end(); ++p) {
		    (*p)->deterministicSample(_chain);
		}
	    }
	}

	double SumMethod::value() const
	{
	    return _x[_i];
	}

	void SumMethod::getLimits(double *lower, double *upper) const
	{
	    vector<StochasticNode *> const &n = _gv->nodes();
	    double li, ui, lj, uj;
	    n[_i]->support(&li, &ui, 1U, _chain);
	    n[_j]->support(&lj, &uj, 1U, _chain);
	    double sum_ij = n[_i]->value(_chain)[0] + n[_j]->value(_chain)[0];
	    *lower = max(li, sum_ij - uj);
	    *upper = min(ui, sum_ij - lj);
	}

	bool SumMethod::isAdaptive() const
	{
	    return true;
	}
	
	void SumMethod::adaptOff()
	{
	    _adapt=false;
	}

	double SumMethod::logDensity() const
	{
	    if (_fast) {
		//Log density depends only on the prior of the two nodes
		//that are currently active
		return _gv->nodes()[_i]->logDensity(_chain, PDF_PRIOR) +
		    _gv->nodes()[_j]->logDensity(_chain, PDF_PRIOR);
	    }
	    else {
		return _gv->logFullConditional(_chain);
	    }
	}

	bool SumMethod::checkAdaptation() const
	{
	    return true;
	}

    } // namespace bugs
} //namespace jags

