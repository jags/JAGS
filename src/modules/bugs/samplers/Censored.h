#ifndef CENSORED_H_
#define CENSORED_H_

#include "ConjugateMethod.h"

#include <graph/ScalarStochasticNode.h>

namespace jags {

    class Graph;

    namespace bugs {
	
	/**
	 * Update Method for censored real-valued distributions.
	 */
	class Censored : public ConjugateMethod {
	    ScalarStochasticNode *_snode;
	public:
	    Censored(SingletonGraphView const *gv);
	    static bool canSample(StochasticNode *snode, Graph const &graph);
	    void update(unsigned int chain, RNG * rng) const;
	};
    
    }
}

#endif /* CENSORED_H_ */
