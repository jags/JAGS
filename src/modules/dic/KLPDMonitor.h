#ifndef KL_PD_MONITOR_H_
#define KL_PD_MONITOR_H_

#include "PDMonitor.h"
#include "KL.h"

class StochasticNode;
class RNG;

namespace dic {

    class KLPDMonitor : public PDMonitor {
	KL const *_kl;
    public:
	KLPDMonitor(StochasticNode const *snode,
		    unsigned int start,  unsigned int thin, KL const *kl);
	void doUpdate();
    };

}

#endif /* KL_PD_MONITOR_H_ */
