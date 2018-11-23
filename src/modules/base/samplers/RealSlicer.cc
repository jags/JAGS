#include <config.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>

#include "RealSlicer.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace base {

    RealSlicer::RealSlicer(SingletonGraphView const *gv, unsigned int chain,
			   double width, unsigned int maxwidth)
	: Slicer(width, maxwidth), _gv(gv), _chain(chain)
    {
	if (!canSample(gv->node()))
	{
	    throwLogicError("Invalid RealSlicer");
	}
	gv->checkFinite(chain);
    }

    bool 
    RealSlicer::canSample(StochasticNode const *node)
    {
	if (node->isDiscreteValued() || node->length() != 1)
	    return false;

	if (!node->fullRank())
	    return false; 

	return true;
    }

    double RealSlicer::value() const
    {
	return _gv->node()->value(_chain)[0];
    }
 
    void RealSlicer::setValue(double value)
    {
	_gv->setValue(&value, 1, _chain);
    }

    void RealSlicer::getLimits(double *lower, double *upper) const
    {
	_gv->node()->support(lower, upper, 1, _chain);
    }

    void RealSlicer::update(RNG *rng)
    {
	if (!updateStep(rng)) {
	    switch(state()) {
	    case SLICER_POSINF:
		throwNodeError(_gv->node(),
			       "Slicer stuck at value with infinite density");
		break;
	    case SLICER_NEGINF:
		throwNodeError(_gv->node(),
			       "Current value is inconsistent with data");
		break;
	    case SLICER_OK:
		break;
	    }
	}
    }

    double RealSlicer::logDensity() const
    {
	return _gv->logFullConditional(_chain);
    }

}}
