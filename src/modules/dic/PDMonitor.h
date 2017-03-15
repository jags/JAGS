#ifndef PD_MONITOR_H_
#define PD_MONITOR_H_

#include <model/Monitor.h>

#include <vector>

namespace jags {

    class StochasticNode;
    struct RNG;
    
namespace dic {

    class PDMonitor : public Monitor {
	std::vector<StochasticNode const *> _snodes;
	std::vector<RNG *> _rngs;
	unsigned int _nrep;
	std::vector<double> _values;
	std::vector<double> _weights;
	double _scale;
	unsigned int _nchain;
    public:
	PDMonitor(std::vector<StochasticNode const *> const &snodes,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep, double scale=1);
	~PDMonitor();
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	bool poolChains() const;
	bool poolIterations() const;
	void update();
	virtual double weight(StochasticNode const *snode,
			      unsigned int ch) const;
    };

}}

#endif /* PD_MONITOR_H_ */
