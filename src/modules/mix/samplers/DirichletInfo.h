#ifndef DIRICHLET_INFO_H_
#define DIRICHLET_INFO_H_

namespace jags {

class StochasticNode;

namespace mix {

    /**
     * Helper class to store information about nodes with a Dirichlet
     * distribution. This is used to invoke special rules for the
     * NormMix sample method.
     */
    struct DirichletInfo 
    {
	unsigned long start;
	unsigned long end;
	unsigned long length;
	double sum;
	double shape;
	DirichletInfo(StochasticNode const *snode, unsigned long start,
		      unsigned int chain);
	double gammaPenalty() const;
    };

}}

#endif /* DIRICHLET_INFO_H_ */
