#ifndef DENSITY_TRACE_H_
#define DENSITY_TRACE_H_

#include <model/Monitor.h>
#include <graph/Node.h>

#include <vector>

namespace jags {
    namespace dic {
	
   	/**
   	 * @short Stores values of density/log density/deviance corresponding to sampled values of a given Node
	 *
	 * Note that this class is used by both NodeDensityMonitorFactory and ObsStochDensMonitorFactory
   	 */
   	class DensityTrace : public Monitor {
 	  protected:
   	    std::vector<Node const *> _nodes;
   	    std::vector<std::vector<double> > _values; // density/log density/deviance corresponding to sampled values
		DensityType const _density_type;  // enum is defined in model/Monitor.h
		std::vector<unsigned int> _dim;
		unsigned int const _nchain;
   	  public:
   	    DensityTrace(std::vector<Node const *> const &nodes, std::vector<unsigned int> dim, 
				DensityType const density_type, std::string const &monitor_name);
   	    void update();
   	    std::vector<double> const &value(unsigned int chain) const;
   	    std::vector<unsigned int> dim() const;
   	    bool poolChains() const;
   	    bool poolIterations() const;
   	};
	
}
}

#endif /* DENSITY_TRACE_H_ */
