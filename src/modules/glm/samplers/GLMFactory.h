#ifndef GLM_FACTORY_H_
#define GLM_FACTORY_H_

#include "GLMMethod.h"
#include <sampler/SamplerFactory.h>

namespace jags {

    class LinkNode;
    class SingletonGraphView;

namespace glm {

    class GLMSampler;
    
    /**
     * @short Abstract factory for GLM samplers 
     * 
     * All factory objects for samplers that handle generalized linear
     * models (GLMs) must recognize the same basic structures in the
     * graph. They differ only by the allowed link function(s) and
     * outcome distribution(s) for the GLM, and the update methods
     * that they use for sampling.
     */
    class GLMFactory : public SamplerFactory
    {
	std::string _name;
	SingletonGraphView * makeView(StochasticNode *snode, 
				      Graph const &graph, bool gibbs) const;
	bool checkDescendants(GraphView const *view) const;
    public:
	GLMFactory(std::string const &name);
	virtual ~GLMFactory();
	/**
	 * Creates a newly-allocated Sampler for the set of candidate nodes,
	 * or a NULL pointer. Sub-classes of GLMFactory only have to
	 * implement the abstract member function newMethod.
	 */
	GLMSampler * makeSampler(std::list<StochasticNode*> const &free_nodes, 
				 Graph const &graph, bool gibbs) const;
	/**
	 * Wraps GLMFactory#makeSampler and returns a single
	 * newly-allocated sampler in a vector.
	 */
	std::vector<Sampler*> 
	    makeSamplers(std::list<StochasticNode*> const &free_nodes, 
			 Graph const &graph) const;
	/**
	 * Checks that an outcome variable in a GLM has the correct
	 * family and link function.
	 *
	 * @param snode Stochastic node representing an outcome variable
	 */
	virtual bool checkOutcome(StochasticNode const *snode) const = 0;
	/**
	 * Function for returning a newly allocated GLMMethod.  This
	 * is called by GLMFactory#makeSampler
	 *
	 * @param view Pointer to a GraphView object for all sampled nodes.
	 *
	 * @param sub_views Vector of pointers to SingletonGraphView
	 * objects with length equal to the number of sampled
	 * nodes. Each sub-view corresponds to a single sampled node.
	 *
	 * @param chain Number of the chain (starting from 0) to which
	 * the sampling method will be applied.
	 *
	 * @param gibbs Flag indicating whether sampling method should do
	 * Gibbs updates, in which nodes are updated individually. If false
	 * then all nodes should be updated in a block.
	 */
	virtual GLMMethod *
	    newMethod(GraphView const *view,
		      std::vector<SingletonGraphView const *> const &sub_views, 
		      unsigned int chain, bool gibbs) const = 0;
	/**
	 * Returns the string provided as the parameter "name" in the
	 * constructor.
	 */
	std::string name() const;
	/**
	 * Returns true if the factory can only work with outcome
	 * variables that are observed. The default implementation
	 * returns false.
	 *
	 * If a sampling method uses auxiliary variables, then
	 * generally it cannot work with unobserved outcomes.  With
	 * auxiliary variable sampling the observed outcome is
	 * considered to be a deterministic function of the auxiliary
	 * variables. When an unobserved outcome variable is updated
	 * (by a sampling method independent of this one) it may be
	 * inconsistent with the current values of the auxiliary
	 * variables, leaving them in an invalid state. An exception
	 * is when the auxiliary variables can be sampled directly,
	 * without reference to the current values.
	 */
	virtual bool fixedOutcome() const;
	/**
	 * Returns true if the factory can only work on a GLM with a
	 * fixed design matrix (i.e. one that does not change between
	 * iterations) and false otherwise.  GLMs with variable design
	 * matrices occur in the context of missing data and
	 * measurement error models.  The default implementation
	 * returns false.
	 *
	 * Note that the same restrictions as for the fixedOutome
	 * function apply for auxiliary variables.
	 */
	virtual bool fixedDesign() const;
	/**
	 * Makes random effects variance samplers.
	 *
	 * Samplers for the variances of random effects in a GLM may
	 * exploit the same data augmentation used by the GLMSampler.
	 */
	void makeRESamplers(std::list<StochasticNode*> const &free_nodes,
			    GLMSampler const *s, Graph const &graph,
			    std::vector<Sampler*> &samplers) const;
    };

}}

#endif /* GLM_FACTORY_H_ */
