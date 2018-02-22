#ifndef DENSITY_TRACE_H_
#define DENSITY_TRACE_H_

#include <model/Monitor.h>
#include <model/NodeArraySubset.h>

#include <vector>

namespace jags {
    namespace dic {

   	/**
   	 * @short Stores values of density/log density/deviance corresponding to sampled values of a given Node
   	 */
   	class DensityTrace : public Monitor {
 	  protected:
   	    NodeArraySubset _subset;
   	    std::vector<std::vector<double> > _values; // density/log density/deviance corresponding to sampled values
		DensityType const _density_type;
   	  public:
   	    DensityTrace(NodeArraySubset const &subset, DensityType const density_type, std::string const &monitor_name);
   	    void update();
   	    std::vector<double> const &value(unsigned int chain) const;
   	    std::vector<unsigned int> dim() const;
   	    bool poolChains() const;
   	    bool poolIterations() const;
   	};
	
}
}

#endif /* DENSITY_TRACE_H_ */
