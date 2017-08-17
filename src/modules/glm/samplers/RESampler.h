#ifndef RE_SAMPLER_H_
#define RE_SAMPLER_H_

#include <sampler/Sampler.h>

namespace jags {

    struct RNG;
    class SingletonGraphView;

    namespace glm {

	class REMethod;

	/**
	 * @short Multi-level sampler for random effects and their precision 
	 *
	 * A vector of REMethod objects is used to update each chain
	 * in parallel.
	 */
	class RESampler : public Sampler
	{
	    SingletonGraphView *_tau;
	    GraphView *_eps;
	    std::vector<SingletonGraphView*> _sub_eps;
	    std::vector<REMethod*> _methods;
	    const std::string _name;
	  public:
	    /**
	     * Constructor.
	     *
	     * @param view Multilevel view of the sample graph, passed
	     * directly to the parent class Sampler, which takes
	     * ownership of it
	     *
	     * @param tau view of the precision parameter and its
	     * descendants. The RESampler takes ownership and deletes when
	     * its destructor is called.
	     *
	     * @param eps view of the random effects and their descendants.
	     * The RESampler takes ownership.
	     *
	     * @param sub_eps vector of singleton views: one for each
	     * random effect. The RESampler takes ownership.
	     *
	     * @param methods Vector of pointers to REMethod
	     * objects, of length equal to the number of chains.  These
	     * must be dynamically allocated, as the
	     * RESampler will take ownership of them, and
	     * will delete them when its destructor is called
	     *
	     * @param name The name of the sampler, which will be returned
	     * by the member function name.
	     */
	    RESampler(GraphView *view,
		      SingletonGraphView *tau, GraphView *eps,
		      std::vector<SingletonGraphView*> sub_eps,
		      std::vector<REMethod*> const &methods,
		      std::string const &name);
	    ~RESampler();
	    void update(std::vector<RNG*> const &rngs);
	    bool isAdaptive() const;
	    void adaptOff();
	    bool checkAdaptation() const;
	    /**
	     * Returns the name of the sampler, as given to the constructor
	     */
	    std::string name() const;
	};

    } // namespace glm

} // namespace jags

#endif /* RE_SAMPLER_H_ */
