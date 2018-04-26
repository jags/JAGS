#ifndef DEVIANCE_TRACE_H_
#define DEVIANCE_TRACE_H_

#include <model/Monitor.h>

namespace jags {

class StochasticNode;

namespace dic {

    class DevianceTrace : public Monitor {
	std::vector<std::vector<double> >  _values; // sampled values
	std::vector<StochasticNode const *> _snodes;
    public:
	DevianceTrace(std::vector<StochasticNode const *> const &nodes);
	std::vector<unsigned long> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void update();
	bool poolChains() const;
	bool poolIterations() const;
    };

}}

#endif /* DEVIANCE_TRACE_H_ */
