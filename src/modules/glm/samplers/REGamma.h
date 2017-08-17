#ifndef RE_GAMMA_H_
#define RE_GAMMA_H_

#include "REMethod.h"
#include "REGammaSlicer.h"

namespace jags {
    namespace glm {

	/**
	 * @short Random effects sampler for gamma precision
	 */
	class REGamma : public REMethod {
	    REGammaSlicer _slicer;
	  public:
	    REGamma(SingletonGraphView const *tau,
		    GraphView const *eps, 
		    std::vector<SingletonGraphView const *> const &veps,
		    std::vector<Outcome *> const &outcomes,
		    unsigned int chain);
	    void updateTau(RNG *rng);
	    void updateSigma(RNG *rng);
	    bool isAdaptive() const;
	    void adaptOff();
	    bool checkAdaptation() const;
	};

    }
}

#endif /* RE_GAMMA_H_ */
