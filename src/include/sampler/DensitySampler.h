#ifndef DENSITY_SAMPLER_H_
#define DENSITY_SAMPLER_H_

#include <sampler/Sampler.h>

class RNG;
class DensityMethod;

/**
 * @short Samples multiple chains in parallel 
 *
 * A DensitySampler uses a vector of DensityMethod objects to update
 * each chain independently.
 */
class DensitySampler : public Sampler
{
    std::vector<DensityMethod*> _methods;
public:
    /**
     * Constructor.
     *
     * @param methods Vector of pointers to  DensityMethod objects,
     * These must be dynamically allocated, as the DensitySampler
     * will take ownership of them, and will delete them when its
     * destructor is called
     */
    DensitySampler(Updater *updater,
		   std::vector<DensityMethod*> const &methods);
    ~DensitySampler();
    void update(std::vector<RNG*> const &rngs);
    bool adaptOff();
    bool isAdaptive() const;
    std::string name() const;
};

#endif /* DENSITY_SAMPLER_H_ */
