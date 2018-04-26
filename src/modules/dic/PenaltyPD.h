#ifndef PENALTY_PD_H_
#define PENALTY_PD_H_

#include <model/Monitor.h>
#include <graph/Node.h>
#include <rng/RNG.h>

#include <vector>


namespace jags {
namespace dic {

    class PenaltyPD : public Monitor {
	protected:
	std::vector<Node const *> const _nodes;
	std::vector<RNG *> _rngs;
	unsigned int _nrep;
	std::vector<double> _values;
	std::vector<unsigned long> const _dim;
	double _scale_cst;
	unsigned long _nchain;
	unsigned int _n;

	// Protected constructor for use by PenaltyPOPT:
	PenaltyPD(std::vector<Node const *> const &nodes,
		  std::vector<unsigned long> const &dim, 
		  std::string const &monitor_name,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep, double scale);
	
    public:
	PenaltyPD(std::vector<Node const *> const &nodes,
		  std::vector<unsigned long> const &dim, 
		  std::string const &monitor_name,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep);

	~PenaltyPD();
	std::vector<unsigned long> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	bool poolChains() const;
	bool poolIterations() const;
	void update();
	};

}}

#endif /* PENALTY_PD_H_ */
