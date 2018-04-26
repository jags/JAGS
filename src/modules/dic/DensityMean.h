#ifndef DENSITY_MEAN_H_
#define DENSITY_MEAN_H_

#include <model/Monitor.h>
#include <graph/Node.h>

#include <vector>

#include "DensityEnums.h"

namespace jags {
    namespace dic {
	
   	/**
   	 * @short Stores running mean values of density/log density/deviance for a given Node
	 *
	 * Note that this class is used by both NodeDensityMonitorFactory and ObsStochDensMonitorFactory
   	 */
   	class DensityMean : public Monitor {
 	  protected:
   	    std::vector<Node const *> const _nodes;
   	    std::vector<std::vector<double> > _values; // density/log density/deviance corresponding to sampled values
		DensityType const _density_type;  // enum is defined in model/Monitor.h
		std::vector<unsigned long> const _dim;
		unsigned int const _nchain;
		unsigned int _n;
   	  public:
   	    DensityMean(std::vector<Node const *> const &nodes, std::vector<unsigned long> const &dim, 
				DensityType const density_type, std::string const &monitor_name);
   	    void update();
   	    std::vector<double> const &value(unsigned int chain) const;
   	    std::vector<unsigned long> dim() const;
   	    bool poolChains() const;
   	    bool poolIterations() const;
   	};
	
}
}

#endif /* DENSITY_MEAN_H_ */
