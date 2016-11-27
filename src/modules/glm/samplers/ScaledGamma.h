#ifndef SCALED_GAMMA_H_
#define SCALED_GAMMA_H_

#include <sampler/SampleMethodNoAdapt.h>

#include <vector>

namespace jags {
    
    class Graph;
    class SingletonGraphView;
    class StochasticNode;
    class RNG;

    namespace glm {
	
	class ScaledGamma : public SampleMethodNoAdapt {
	    SingletonGraphView const *_gv;
	    unsigned int _chain;
	    std::vector<double> _coef;
	    double _a;
	    bool _fast;
	    void calCoef();
	  public:
	    ScaledGamma(SingletonGraphView const *gv, unsigned int chain);
	    static bool canSample(StochasticNode *snode, Graph const &graph);
	    void update(RNG *rng);
	};
	
    }
}

#endif /* SCALED_GAMMA_H_ */
