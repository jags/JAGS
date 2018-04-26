#ifndef VECTOR_STOCHASTIC_NODE_H_
#define VECTOR_STOCHASTIC_NODE_H_

#include <graph/StochasticNode.h>

namespace jags {

class VectorDist;

/**
 * @short Vector-valued Node defined by the BUGS-language operator ~
 */
class VectorStochasticNode : public StochasticNode {
    VectorDist const * const _dist;
    std::vector<unsigned long> _lengths;
    void sp(double *lower, double *upper, unsigned long length,
	    unsigned int chain) const;
public:
    /**
     * Constructs a new StochasticNode given a vector distribution and
     * a vector of parent nodes, considered as parameters to the
     * distribution.
     */
    VectorStochasticNode(VectorDist const *dist, unsigned int nchain,
			 std::vector<Node const *> const &parameters,
			 Node const *lower, Node const *upper);
    double logDensity(unsigned int chain, PDFType type) const;
    void randomSample(RNG *rng, unsigned int chain);
    void truncatedSample(RNG *rng, unsigned int chain,
			 double const *lower, double const *upper);
    bool checkParentValues(unsigned int chain) const;
    //StochasticNode *clone(std::vector<Node const *> const &parents,
    //Node const *lower, Node const *upper) const;
    unsigned long df() const;
    double KL(unsigned int chain1, unsigned int chain2, RNG *rng,
	      unsigned int nrep) const;
};

} /* namespace jags */

#endif /* VECTOR_STOCHASTIC_NODE_H_ */

