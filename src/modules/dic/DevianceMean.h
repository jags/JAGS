#ifndef DEVIANCE_MEAN_H_
#define DEVIANCE_MEAN_H_

#include <model/Monitor.h>

namespace jags {

class StochasticNode;

namespace dic {

    class DevianceMean : public Monitor {
	std::vector<double>  _values; 
	std::vector<StochasticNode const *> _snodes;
	unsigned int _n;
    public:
	DevianceMean(std::vector<StochasticNode const *> const &nodes);
	std::vector<unsigned long> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void update();
	bool poolChains() const;
	bool poolIterations() const;
    };

}}

#endif /* DEVIANCE_MEAN_H_ */
