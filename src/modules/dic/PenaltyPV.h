#ifndef PENALTY_PV_H_
#define PENALTY_PV_H_

#include <model/Monitor.h>
#include <graph/Node.h>

#include <vector>

namespace jags {
namespace dic {

   	class PenaltyPV : public Monitor {
   	    std::vector<Node const *> _nodes;
		double _mean;
		double _mm;
		std::vector<double> _pv;
		std::vector<unsigned int> _dim;
		unsigned int const _nchain;
		unsigned int _n;
   	  public:
   	    PenaltyPV(std::vector<Node const *> const &nodes,
		      std::string const &monitor_name);
   	    void update();
   	    std::vector<double> const &value(unsigned int chain) const;
   	    std::vector<unsigned int> dim() const;
   	    bool poolChains() const;
   	    bool poolIterations() const;
   	};

}}

#endif /* PENALTY_PV_H_ */
