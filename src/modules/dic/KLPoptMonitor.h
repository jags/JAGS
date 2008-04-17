#ifndef KL_POPT_MONITOR_H_
#define KL_POPT_MONITOR_H_

#include "PoptMonitor.h"
#include "KL.h"

class StochasticNode;

namespace dic {

    class KLPoptMonitor : public PoptMonitor {
	KL const *_kl;
    public:
	KLPoptMonitor(StochasticNode const *snode,
		      unsigned int start,  unsigned int thin, KL const *kl);
	void doUpdate();
    };

}

#endif /* KL_POPT_MONITOR_H_ */
