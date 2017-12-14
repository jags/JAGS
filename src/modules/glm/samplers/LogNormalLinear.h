#ifndef LOG_NORMAL_LINEAR_H_
#define LOG_NORMAL_LINEAR_H_

#include <config.h>
#include "Outcome.h"

#include <graph/StochasticNode.h>

namespace jags {

    class StochasticNode;

    namespace glm {
	
	class LogNormalLinear : public Outcome
	{
	    double const &_value;
	    double const &_precision;
	public:
	    LogNormalLinear(StochasticNode const *snode, unsigned int chain);
	    double value() const;
	    double precision() const;
	    static bool canRepresent(StochasticNode const *snode);
	    bool fixedb() const;
	    bool fixedA() const;
	};

    }
}

#endif /* LOG_NORMAL_LINEAR_H_ */
