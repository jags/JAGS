#ifndef GLM_SAMPLER_H_
#define GLM_SAMPLER_H_

#include <sampler/Sampler.h>

#include <vector>

namespace jags {

    class GraphView;
    class SingletonGraphView;
    struct RNG;
    
namespace glm {

    class GLMMethod;
    class REFactory2;
    
    /**
     * @short Base class for GLM samplers.
     */
    class GLMSampler : public Sampler
    {
	GraphView const *_view;
	std::vector<SingletonGraphView*> _sub_views;
	std::vector<GLMMethod*> _methods;
	std::string _name;
	friend REFactory2;
    public:
	/**
	 * Constructor.
	 *
	 * @param view Pointer to a GraphView object for all sampled nodes.
	 *
	 * @param sub_views Vector of pointers to SingletonGraphView
	 * objects with length equal to the number of sampled
	 * nodes. Each sub-view corresponds to a single sampled node.
	 * The GLMSampler object takes ownership of these sub-views
	 * and deletes them when its destructor is called.
	 *
	 * @param methods Vector of sampling methods
	 */
	GLMSampler(GraphView *view, 
		   std::vector<SingletonGraphView*> const &sub_views,
		   std::vector<GLMMethod*> const &methods,
		   std::string const &name);
	/**
	 * Destructor
	 * 
	 * Deletes the sub-views and the methods passed to the constructor.
	 */
	~GLMSampler();
	void update(std::vector<RNG*> const &rngs);
	bool isAdaptive() const;
	void adaptOff();
	bool checkAdaptation() const;
	std::string name() const;
	/*
	  Gives access to the vector of GLMMethod objects used by the
	  sampler.
	 */
	std::vector<GLMMethod*> const &methods();

    }; 

}}

#endif /* GLM_SAMPLER_H_ */
