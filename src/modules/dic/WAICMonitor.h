#ifndef WAIC_MONITOR_H_
#define WAIC_MONITOR_H_

#include <model/Monitor.h>

#include <vector>

namespace jags {

    class StochasticNode;
    struct RNG;
    
    namespace dic {

	class WAICMonitor : public Monitor {
	    std::vector<StochasticNode const *> _snodes;
	    unsigned int _nchain;
	    std::vector<std::vector<double> > _mlik;
	    std::vector<std::vector<double> > _vlik;
	    unsigned int _n;

	public:
	    WAICMonitor(std::vector<StochasticNode const *> const &snodes);
	    ~WAICMonitor();
	    std::vector<unsigned int> dim() const;
	    std::vector<double> const &value(unsigned int chain) const;
	    bool poolChains() const;
	    bool poolIterations() const;
	    void update();
	};

    }
}

#endif /* WAIC_MONITOR_H_ */
