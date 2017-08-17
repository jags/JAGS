#ifndef M_NORMAL_LINEAR_H_
#define M_NORMAL_LINEAR_H_

#include <config.h>
#include "Outcome.h"

#include <graph/StochasticNode.h>

namespace jags {

class StochasticNode;

namespace glm {

    class MNormalLinear : public Outcome
    {
	double const *_value;
	double const *_mean;
	double const *_precision;
      public:
	MNormalLinear(StochasticNode const *snode, unsigned int chain);
	double value() const;
	double precision() const;
	double const *vvalue() const;
	double const *vmean() const;
	double const *vprecision() const;
	static bool canRepresent(StochasticNode const *snode);
    };

}}

#endif /* M_NORMAL_LINEAR_H_ */
