#ifndef KL_POPT_MONITOR_H_
#define KL_POPT_MONITOR_H_

#include <model/Monitor.h>
#include <graph/StochasticNode.h>

#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class KLPoptMonitor : public Monitor {
	KL const *_kl;
	std::vector<std::vector<double const *> > _par;
	std::vector<double> _weights; // weights for each chain
	std::vector<double> _values; // sampled values
    public:
	PoptMonitor(StochasticNode const *snode,
		    unsigned int start,  unsigned int thin, 
		    std::vector<RNG *> const &rngs, unsigned int nrep); 
	unsigned int nchain() const;
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void doUpdate();
	void reserve(unsigned int niter);
	SArray dump() const;
    };

}

#endif /* KL_POPT_MONITOR_H_ */
