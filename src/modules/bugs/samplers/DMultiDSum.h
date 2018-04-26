#ifndef DMULTI_DSUM_H_
#define DMULTI_DSUM_H_

#include "RWDSum.h"

class StochasticNode;

namespace jags {
namespace bugs {

/**
 * @short Sample parents of dsum nodes
 */
class DMultiDSum : public RWDSum
{
public:
    DMultiDSum(GraphView const *gv, unsigned int chain);
    void step(std::vector<double> &value, unsigned int long,
	      unsigned long ncol, double s, RNG *rng) const;
};

}}

#endif /* DMULTI_DSUM_H_ */
