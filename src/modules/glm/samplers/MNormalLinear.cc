#include <config.h>

#include "MNormalLinear.h"
#include "Classify.h"

#include <graph/StochasticNode.h>

namespace jags {
    namespace glm {

	MNormalLinear::MNormalLinear(StochasticNode const *snode,
				     unsigned int chain)
	    : Outcome(snode, chain),
	      _value(snode->value(chain)),
	      _mean(snode->parents()[0]->value(chain)),
	      _precision(snode->parents()[1]->value(chain))
	{
	}
	
	double MNormalLinear::value() const 
	{
	    return 0;
	}

	double MNormalLinear::precision() const 
	{
	    return 0;
	}

	bool MNormalLinear::canRepresent(StochasticNode const *snode)
	{
	    return snode->distribution()->name() == "dmnorm" &&
		getLink(snode) == LNK_LINEAR;
	}

	double const *MNormalLinear::vvalue() const
	{
	    return _value;
	}

	double const *MNormalLinear::vmean() const
	{
	    return _mean;
	}

	double const *MNormalLinear::vprecision() const
	{
	    return _precision;
	}

     
	
    }}
