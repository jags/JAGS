#include <config.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>

#include "BinomSlicer.h"

#include <cmath>

using std::vector;
using std::string;
using std::log;

namespace jags {
    namespace bugs {

    BinomSlicer::BinomSlicer(SingletonGraphView const *gv, unsigned int chain,
			   double width, unsigned int maxwidth)
	: Slicer(width, maxwidth), _gv(gv), _chain(chain)
    {
	gv->checkFinite(chain);
    }

    bool 
    BinomSlicer::canSample(StochasticNode *node, Graph const &graph)
    {
	if (node->isDiscreteValued() || node->length() != 1)
	    return false;

	if (node->df() == 0)
	    return false; 

	SingletonGraphView gv(node, graph);
	vector<StochasticNode *> const &schild = gv.stochasticChildren();
	for (unsigned int i = 0; i < schild.size(); ++i) {
	    if (schild[i]->distribution()->name() != "dbin") {
		return false;
	    }
	    if (!schild[i]->parents()[1]->isFixed()) {
		return false;
	    }
	}

	return true;
    }

    double BinomSlicer::value() const
    {
	return _gv->node()->value(_chain)[0];
    }
 
    void BinomSlicer::setValue(double value)
    {
	_gv->setValue(&value, 1, _chain);
    }

    void BinomSlicer::getLimits(double *lower, double *upper) const
    {
	_gv->node()->support(lower, upper, 1, _chain);
    }

    void BinomSlicer::update(RNG *rng)
    {
	if (!updateStep(rng)) {
	    switch(state()) {
	    case SLICER_POSINF:
		throwNodeError(_gv->node(),
			       "Slicer stuck at value with infinite density");
	    case SLICER_NEGINF:
		throwNodeError(_gv->node(),
			       "Current value is inconsistent with data");
	    case SLICER_OK:
		break;
	    }
	}
    }

    double BinomSlicer::logDensity() const
    {
	//Override expensive log likelihood calculations for binomial
	double loglik =  _gv->logPrior(_chain);
	vector<StochasticNode *> const &schild = _gv->stochasticChildren();
	for (unsigned int i = 0; i < schild.size(); ++i) {
	    double y = schild[i]->value(_chain)[0];
	    double p = schild[i]->parents()[0]->value(_chain)[0];
	    double n = schild[i]->parents()[1]->value(_chain)[0];
	    if (y==0) {
		//FIXME: Need C++11 for log1p
		loglik += n*log(1-p);
	    }
	    else if (y==n) {
		loglik += y*log(p);
	    }
	    else {
		loglik += y*log(p)+(n-y)*log(1-p);
	    }
	}
        return loglik;
    }

    }
}
