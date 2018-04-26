#ifndef IMMUTABLE_SAMPLER_H_
#define IMMUTABLE_SAMPLER_H_

#include <sampler/Sampler.h>

namespace jags {

    struct RNG;
    class ImmutableSampleMethod;

    /**
     * @short Samples multiple chains in parallel
     *
     * Uses an ImmutableSampleMethod object to update each chain
     * independently.
     *
     * In contrast to the MutableSampler, only one sample
     * method is required to update all chains.
     */
    class ImmutableSampler : public Sampler
    {
	ImmutableSampleMethod const * const _method;
	const unsigned int _nchain;
	const std::string _name;
      public:
	/**
	 * Constructor.
	 *
	 * @param gv View of the sample graph, passed directly to the
	 * parent class Sampler.
	 *
	 * @param method Pointer to an ImmutableSampleMethod. This
	 * must be dynamically allocated, as the
	 * ImmutableSampler will take ownership of it, and
	 * will delete it when the destructor is called.
	 *
	 * @param name Name of the sampler, which will be returned
	 * by the member function name.
	 */
	ImmutableSampler(GraphView *gv, ImmutableSampleMethod* method,
			 std::string const &name);
	~ImmutableSampler();
	/**
	 * Updates all chains in parallel using the update method
	 * provided to the constructor
	 */
	void update(std::vector<RNG*> const &rngs);
	/**
	 * The sampler is not adaptive
	 */
	bool isAdaptive() const;
	/**
	 * This does nothing
	 */
	void adaptOff();
	/**
	 * This always returns true
	 */
	bool checkAdaptation() const;
	/**
	 * Returns the name of the sampler, as given to the constructor.b
	 */
	std::string name() const;
    };

} /* namespace jags */

#endif /* IMMUTABLE_SAMPLER_H_ */
