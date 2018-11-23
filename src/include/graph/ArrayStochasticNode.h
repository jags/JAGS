#ifndef ARRAY_STOCHASTIC_NODE_H_
#define ARRAY_STOCHASTIC_NODE_H_

#include <graph/StochasticNode.h>

namespace jags {

class ArrayDist;

/**
 * @short Array-valued Node defined by the BUGS-language operator ~
 */
class ArrayStochasticNode : public StochasticNode {
    ArrayDist const * const _dist;
    std::vector<std::vector<unsigned long> > _dims;
    void sp(double *lower, double *upper, unsigned long length, 
	    unsigned int chain) const;
public:
    /**
     * Construct
     * distribution and an array of parent nodes, considered as
     * parameters to the distribution.
     */
    ArrayStochasticNode(ArrayDist const *dist, unsigned int nchain,
			std::vector<Node const *> const &parameters);
    double logDensity(unsigned int chain, PDFType type) const;
    void randomSample(RNG *rng, unsigned int chain);
    bool checkParentValues(unsigned int chain) const;
    //StochasticNode *clone(std::vector<Node const *> const &parents,
    //Node const *lower, Node const *upper) const;
    double KL(unsigned int chain1, unsigned int chain2, RNG *rng, unsigned int nrep) const;
};

} /* namespace jags */

#endif /* ARRAY_STOCHASTIC_NODE_H_ */

