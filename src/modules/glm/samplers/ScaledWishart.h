#ifndef SCALED_WISHART_H_
#define SCALED_WISHART_H_

#include <sampler/SampleMethodNoAdapt.h>

#include <vector>

namespace jags {

    class Graph;
    class SingletonGraphView;
    class StochasticNode;
    class RNG;
    
    namespace glm {

	class ScaledWishart : public SampleMethodNoAdapt {
	    SingletonGraphView const *_gv;
	    unsigned int _chain;
	    std::vector<double> _a;
	  public:
	    ScaledWishart(SingletonGraphView const *gv, unsigned int chain);
	    void update(RNG *rng);
	    static bool canSample(StochasticNode *snode, Graph const &graph);
	};

    }
}

#endif /* SCALED_WISHART_H_ */
