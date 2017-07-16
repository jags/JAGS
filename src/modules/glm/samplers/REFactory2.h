#ifndef RE_FACTORY2_H_
#define RE_FACTORY2_H_

#include <sampler/Sampler.h>

#include <set>
#include <string>
#include <list>

namespace jags {

    class SingletonGraphView;
    class GraphView;
    class Graph;
    
    namespace glm {

	class REMethod2;
	class GLMSampler;
	class GLMMethod;
	
	class REFactory2 
	{
	    const std::string _name;
	  public:
	    bool checkTau(SingletonGraphView const *tau,
			  GraphView const *glmview) const;

	    REFactory2(std::string const &name);
	    virtual ~REFactory2();
	    Sampler *makeSampler(std::list<StochasticNode*> const &free_nodes,
				 std::set<StochasticNode*> &used_nodes,
				 GLMSampler const *glmsampler,
				 Graph const &graph) const;
	    
	    virtual bool canSample(StochasticNode *snode) const = 0;
	    virtual REMethod2 *
		newMethod(SingletonGraphView const *tau,
			  GLMMethod const *randef) const = 0;
	};
	
    }
}

#endif /* RE_FACTORY_H_ */
