#include <config.h>
#include <graph/StochasticNode.h>
#include <graph/NodeError.h>
#include <distribution/Distribution.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include "DiscreteSlicer.h"

#include <cmath>

using std::floor;
using std::vector;
using std::string;

namespace jags {
namespace base {

    DiscreteSlicer::DiscreteSlicer(SingletonGraphView const *gv, 
				   unsigned int chain, 
				   double width, unsigned int ndoubles)
	: Slicer(width, ndoubles), _gv(gv), _chain(chain), _x(0)
    {
	if (!canSample(gv->node())) {
	    throwLogicError("Invalid DiscreteSlicer");
	}

	//Check finite likelihood and copy current value
	gv->checkFinite(chain);
	_x = _gv->node()->value(chain)[0];
    }

    bool DiscreteSlicer::canSample(StochasticNode const *node)
    {
	if (!node->isDiscreteValued() || node->length() != 1)
	    return false;
	
	return true;
    }

    void DiscreteSlicer::setValue(double x)
    {
	_x = x;
	x = floor(x);
	_gv->setValue(&x, 1, _chain);
    }
  
    double DiscreteSlicer::value() const
    {
	return _x;
    }

    void DiscreteSlicer::getLimits(double *lower, double *upper) const
    {
        _gv->node()->support(lower, upper, 1, _chain);
	*upper += 1;
    }
    
    void DiscreteSlicer::update(RNG *rng)
    {
	if (!updateDouble(rng)) {
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

    double DiscreteSlicer::logDensity() const
    {
	return _gv->logFullConditional(_chain);
    }

}}
