#ifndef RE_FACTORY_H_
#define RE_FACTORY_H_

#include <sampler/SamplerFactory.h>

namespace jags {

    class SingletonGraphView;
    class GraphView;
    
    namespace glm {

	class REMethod;
	class Outcome;
	
	/**
	 * @short Factory for scaled precision parameters 
	 *
	 * @see ScaledGamma
	 */
	class REFactory : public SamplerFactory
	{
	    std::string _name;
	    bool checkOutcome(StochasticNode const *snode) const;
	    bool checkTau(SingletonGraphView const *tau) const;
	    bool checkEps(GraphView const *eps) const;
	  public:
	    REFactory(std::string const &name);
	    ~REFactory();
	    std::vector<Sampler*> 
		makeSamplers(std::list<StochasticNode*> const &nodes, 
			     Graph const &graph) const;
	    Sampler *makeSampler(std::list<StochasticNode*> const &free_nodes, 
				 Graph const &graph) const;
	    std::string name() const;
	    
	    virtual bool canSample(StochasticNode *snode) const = 0;
	    virtual REMethod *
		newMethod(SingletonGraphView const *tau,
			  GraphView const *eps,
			  std::vector<SingletonGraphView const *> const & veps,
			  std::vector<Outcome*> const &outcomes,
			  unsigned int chain) const = 0;
		
	};
	
    }
}

#endif /* RE_FACTORY_H_ */
