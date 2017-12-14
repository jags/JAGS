#include <config.h>

#include <cmath>

#include "LogNormalLinear.h"
#include "Classify.h"

#include <graph/StochasticNode.h>

using std::log;

namespace jags {
    namespace glm {
	
	LogNormalLinear::LogNormalLinear(StochasticNode const *snode,
					 unsigned int chain)
	    : Outcome(snode, chain),
	      _value(snode->value(chain)[0]),
	      _precision(snode->parents()[1]->value(chain)[0])
	{
	}
    
	double LogNormalLinear::value() const 
	{
	    return log(_value);
	}

	double LogNormalLinear::precision() const 
	{
	    return _precision;
	}

	bool LogNormalLinear::canRepresent(StochasticNode const *snode)
	{
	    return getFamily(snode) == GLM_LNORMAL &&
		getLink(snode) == LNK_LINEAR;
	}

	bool LogNormalLinear::fixedb() const
	{
	    return true;
	}

	bool LogNormalLinear::fixedA() const
	{
	    return true;
	}
    

    } /* namespace glm */
} /* namespace jags */
