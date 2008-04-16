#ifndef KL_POPT_MONITOR_H_
#define KL_POPT_MONITOR_H_

#include "PoptMonitor.h"

class StochasticNode;

namespace dic {

    class KL;

    class KLPoptMonitor : public PoptMonitor {
	KL const *_kl;
    public:
	KLPoptMonitor(StochasticNode const *snode,
		      unsigned int start,  unsigned int thin, KL const *kl);
	unsigned int nchain() const;
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void doUpdate();
	void reserve(unsigned int niter);
	SArray dump() const;
    };

}

#endif /* KL_POPT_MONITOR_H_ */
