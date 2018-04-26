#ifndef DISCRETE_DSUM_H_
#define DISCRETE_DSUM_H_

#include "RWDSum.h"

class StochasticNode;

namespace jags {
namespace bugs {

/**
 * @short Sample parents of dsum nodes
 */
class DiscreteDSum : public RWDSum
{
public:
    DiscreteDSum(GraphView const *gv, unsigned int chain);
    void step(std::vector<double> &value, unsigned long nrow,
	      unsigned long ncol, double s, RNG *rng) const;
};

}}

#endif /* DISCRETE_DSUM_H_ */
