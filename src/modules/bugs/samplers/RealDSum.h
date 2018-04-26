#ifndef REAL_DSUM_H_
#define REAL_DSUM_H_

#include "RWDSum.h"

#include <vector>

class StochasticNode;
class GraphView;

namespace jags {
namespace bugs {

/**
 * @short Sample parents of dsum nodes
 */
class RealDSum : public RWDSum
{
public:
    RealDSum(GraphView const *gv, unsigned int chain);
    void step(std::vector<double> &value, unsigned long nrow, 
	      unsigned long ncol, double s, RNG *rng) const;
};

}}

#endif /* REAL_DSUM_H_ */
