#ifndef IWLS_FACTORY_H_
#define IWLS_FACTORY_H_

#include "GLMFactory.h"

namespace jags {
namespace glm {

    /**
     * @short Factory object for iteratively weighted least squares
     */
    class IWLSFactory : public GLMFactory
    {
    public:
	IWLSFactory();
	/**
	 * Checks that the outcome variable is a normal, binomial, or
	 * Poisson variable. For normal outcomes, checks that the link
	 * is identity.
	 * 
	 * The IWLS method is potentially open to other distributional
	 * families, but in the BUGS language, these are not
	 * parameterized in terms of their mean and precision.
	 */
	bool checkOutcome(StochasticNode const *snode) const;
	/**
	 * Returns a newly allocated object of class IWLS
	 */
	GLMMethod *newMethod(GraphView const *view, 
			     std::vector<SingletonGraphView const *> const &sub_views,
			     unsigned int chain, bool gibbs) const;
	/**
	 * Returns false if any parents of the candidate node are
	 * unobserved.  The IWLS method relies on an asymptotic
	 * approximation which holds only for fixed effects models.
	 * In order to exclude random effects, we reject candidate
	 * nodes that have unobserved parents.
	 */
	bool canSample(StochasticNode const *snode) const;
	/**
	 * Returns true. The IWLS method relies on an asymptotic
	 * approximation that holds only for fixed effects models.  In
	 * order to exclude random effects, we require the design
	 * matrix to be fixed. This is because it is possible to
	 * reparameterize a random effect as a fixed effect with a
	 * varying coefficient.
	 */
	bool fixedDesign() const;
    };

}}

#endif /* IWLS_FACTORY_H_ */
