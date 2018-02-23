#include "ObsStochDensMonitorFactory.h"
#include "DensityTrace.h"
#include "DensityMean.h"
#include "DensityVariance.h"
#include "DensityTotal.h"
#include "DensityPoolMean.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <graph/StochasticNode.h>
#include <sarray/RangeIterator.h>

#include <set>

using std::set;
using std::string;
using std::vector;

namespace jags {
namespace dic {

    Monitor *ObsStochDensMonitorFactory::getMonitor(string const &name, 
						Range const &range,
						BUGSModel *model,
						string const &type,
						string &msg)
    {
		
		if (name != "deviance")
		    return 0;
		if (!isNULL(range)) {
		    msg = "cannot monitor a subset of the observed stochastic nodes using the name deviance";
		    return 0;
		}
		
		/* Work out the precise type of monitor */
		
		// enums declared in model/Monitor.h:
		MonitorType monitor_type; 
		DensityType density_type; 
		// Used to help resolve aliases:
		string monitor_name = "";
		
		if (type == "density_trace") {
			monitor_type = TRACE;
			density_type = DENSITY;
			monitor_name.assign("density_trace");
		}
		else if (type == "density_mean") {
			monitor_type = MEAN;
			density_type = DENSITY;
			monitor_name.assign("density_mean");
		}
		else if (type == "density_variance") {
			monitor_type = VARIANCE;
			density_type = DENSITY;
			monitor_name.assign("density_variance");
		}
		else if (type == "density_var") {
			monitor_type = VARIANCE;
			density_type = DENSITY;
			monitor_name.assign("density_variance");
		}
		else if (type == "density_total") {
			monitor_type = TOTAL;
			density_type = DENSITY;
			monitor_name.assign("density_total");
		}
		else if (type == "density_poolmean") {
			monitor_type = POOLMEAN;
			density_type = DENSITY;
			monitor_name.assign("density_poolmean");
		}
		else if (type == "logdensity_trace") {
			monitor_type = TRACE;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_trace");
		}
		else if (type == "logdensity_mean") {
			monitor_type = MEAN;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_mean");
		}
		else if (type == "logdensity_variance") {
			monitor_type = VARIANCE;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_variance");
		}
		else if (type == "logdensity_var") {
			monitor_type = VARIANCE;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_variance");
		}
		else if (type == "logdensity_total") {
			monitor_type = TOTAL;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_total");
		}
		else if (type == "logdensity_poolmean") {
			monitor_type = POOLMEAN;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_poolmean");
		}
		else if (type == "deviance_trace") {
			monitor_type = TRACE;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_trace");
		}
		else if (type == "deviance_mean") {
			monitor_type = MEAN;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_mean");
		}
		else if (type == "deviance_variance") {
			monitor_type = VARIANCE;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_variance");
		}
		else if (type == "deviance_var") {
			monitor_type = VARIANCE;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_variance");
		}
		else if (type == "deviance_total") {
			monitor_type = TOTAL;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_total");
		}
		else if (type == "deviance_poolmean") {
			monitor_type = POOLMEAN;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_poolmean");
		}
		else if (type == "trace") {
			// The equivalent of monitor('deviance', type='trace') in DevianceMonitorFactory:
			monitor_type = TOTAL;
			density_type = DEVIANCE;
			monitor_name.assign("trace");  // Note: name is for backwards-compatibility with DevianceMonitorFactory
		}
		else if (type == "mean") {
			// The equivalent of monitor('deviance', type='mean') in DevianceMonitorFactory:
			monitor_type = POOLMEAN;
			density_type = DEVIANCE;
			monitor_name.assign("mean");  // Note: name is for backwards-compatibility with DevianceMonitorFactory
		}
		else {
			// If not listed above:
			return 0;
		}
		
		vector<StochasticNode *> const &snodes = model->stochasticNodes();
		vector<Node const *> observed_snodes;
		for (unsigned int i = 0; i < snodes.size(); ++i) {
		    if (snodes[i]->isFixed()) {
				// Implicit up-cast to Node from StochasticNode:
				observed_snodes.push_back(snodes[i]);
		    }
		}
		if (observed_snodes.empty()) {
		    msg = "There are no observed stochastic nodes";
		    return 0;
		}


		/* Create the correct subtype of monitor */

		Monitor *m = 0;

		if (monitor_type == TRACE) {
			m = new DensityTrace(observed_snodes, density_type, monitor_name);
		}
		else if (monitor_type == MEAN) {
			m = new DensityMean(observed_snodes, density_type, monitor_name);
		}
		else if (monitor_type == VARIANCE) {
			m = new DensityVariance(observed_snodes, density_type, monitor_name);
		}
		else if (monitor_type == TOTAL) {
			m = new DensityTotal(observed_snodes, density_type, monitor_name);
		}
		else if (monitor_type == POOLMEAN) {
			m = new DensityPoolMean(observed_snodes, density_type, monitor_name);
		}
		else {
			throw std::logic_error("Unimplemented MonitorType in ObsStochDensMonitorFactory");
		}
		
		/* Set name attributes */

		m->setName("deviance");

		// TOTAL is the only one summarised between variables
		if (monitor_type == TOTAL) {
		    m->setElementNames(vector<string>(1, monitor_name));
		}
		else {
		    vector<string> onames(observed_snodes.size());
		    for (unsigned int i = 0; i < observed_snodes.size(); ++i) {
				onames[i] = model->symtab().getName(observed_snodes[i]);
		    }
		    m->setElementNames(onames);
		}
		
		return m;
		
    }

    string ObsStochDensMonitorFactory::name() const
    {
		return "dic::ObsStochDensity";
    }
	
}
}
