#ifndef SHIFTED_MULTINOMIAL_H_
#define SHIFTED_MULTINOMIAL_H_

#include "ConjugateMethod.h"

namespace jags {

    class Graph;

namespace bugs {

    class ShiftedMultinomial : public ConjugateMethod {
	std::vector<unsigned long> _index;
    public:
	ShiftedMultinomial(SingletonGraphView const *gv);
	static bool canSample(StochasticNode *snode, Graph const &graph);
	void update(unsigned int chain, RNG *rng) const;
    };

}}

#endif /* SHIFTED_MULTINOMIAL_H_ */
