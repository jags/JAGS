#ifndef CONJUGATE_DIRICHLET_H_
#define CONJUGATE_DIRICHLET_H_

#include "ConjugateMethod.h"

namespace jags {

    class Graph;

namespace bugs {

/**
 * @short Conjugate sampler for Dirichlet distribution.  
 *
 * A conjugate sampler may be created for the stochastic node "snode"
 * if it has a Dirichlet distribution and its marginal stochastic
 * children have either a categorical distribution (dcat) or a
 * multinomial distribution (dmult) with snode as the probability
 * vector.  In the case of a multinomial distribution, the sample size
 * must not depend on snode.  The immediate deterministic children
 * must be either aggregate nodes or mixture nodes.
 */
class ConjugateDirichlet : public ConjugateMethod {
    bool _mix;
    const std::vector<unsigned long> _tree;
    std::vector<std::vector<unsigned long> > _offsets;
    std::vector<unsigned long> _leaves;
    bool isActiveLeaf(unsigned long index, unsigned int chain) const;
    bool isActiveTree(unsigned long index, unsigned int chain) const;
public:
    ConjugateDirichlet(SingletonGraphView const *gv);
    void update(unsigned int chain, RNG *rng) const;
    static bool canSample(StochasticNode *snode, Graph const &graph);
};

}}

#endif /* CONJUGATE_DIRICHLET_H_ */
