#include "DensityEnums.h"
#include "ObsStochDensMonitorFactory.h"
#include "DensityTrace.h"
#include "DensityMean.h"
#include "DensityVariance.h"
#include "DensityTotal.h"
#include "DensityPoolMean.h"
#include "DensityPoolVariance.h"

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
		
		/* Note: the string _observed_stochastic_ is not a valid node
		name so is gauranteed not to clash with a node.
		If changing it then see also scanner.ll and NodeDensityMonitorFactory.cc
		And also parse.varname, waic.samples etc in rjags */
		if (name != "_observed_stochastic_"
			// If replacing DevianceMonitorFactory:
			 && !(name == "deviance" && type == "trace")
			 && !(name == "deviance" && type == "mean"))
		    return 0;
		if (!isNULL(range)) {
		    msg = "cannot monitor a subset of all the observed stochastic nodes - use the specific node names and subsets instead";
		    return 0;
		}
		
		/* Work out the precise type of monitor */
		
		// enums declared in model/Monitor.h:
		MonitorType monitor_type; 
		DensityType density_type; 
		// Used to help resolve aliases:
		string monitor_name = "";
		
		/* The first 2 duplicate DevianceMonitorFactory but the DevianceMonitorFactory
		   takes precedence in dic.cc */
		if (name == "deviance" && type == "trace") {
			// Disable this once it is being used:
			printf("NOTE: Using deviance trace monitor from ObsStochDensMonitorFactory\n");
			// The equivalent of monitor('deviance', type='trace') in DevianceMonitorFactory:
			monitor_type = TOTAL;
			density_type = DEVIANCE;
			// Note: backwards-compatibility with DevianceMonitorFactory
			monitor_name.assign("trace");
		}
		else if (name == "deviance" && type == "mean") {
			// Disable this once it is being used:
			printf("NOTE: Using deviance mean monitor from ObsStochDensMonitorFactory\n");
			// The equivalent of monitor('deviance', type='mean') in DevianceMonitorFactory:
			monitor_type = POOLMEAN;
			density_type = DEVIANCE;
			// Note: backwards-compatibility with DevianceMonitorFactory
			monitor_name.assign("mean");
		}
		else if (type == "density_trace") {
			monitor_type = TRACE;
			density_type = DENSITY;
			monitor_name.assign("density_trace");
		}
		else if (type == "density_mean") {
			monitor_type = MEAN;
			density_type = DENSITY;
			monitor_name.assign("density_mean");
		}
		else if (type == "density_variance" || type == "density_var") {
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
		else if (type == "density_poolvariance" || type == "density_poolvar") {
			monitor_type = POOLVARIANCE;
			density_type = DENSITY;
			monitor_name.assign("density_poolvariance");
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
		else if (type == "logdensity_variance" || type == "logdensity_var") {
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
		else if (type == "logdensity_poolvariance" || type == "logdensity_poolvar") {
			monitor_type = POOLVARIANCE;
			density_type = LOGDENSITY;
			monitor_name.assign("logdensity_poolvariance");
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
		else if (type == "deviance_variance" || type == "deviance_var") {
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
		else if (type == "deviance_poolvariance" || type == "deviance_poolvar") {
			monitor_type = POOLVARIANCE;
			density_type = DEVIANCE;
			monitor_name.assign("deviance_poolvariance");
		}
		else {
			// If not listed above:
			return 0;
		}
		
		vector<Node const *> const &observed_snodes = model->observedStochasticNodes();
		
		if (observed_snodes.empty()) {
		    msg = "There are no observed stochastic nodes";
		    return 0;
		}


		/* Create the correct subtype of monitor */

		Monitor *m = 0;
		
		// There is only ever a single dimension of variables to worry about:
		vector<unsigned int> dim;
		dim.push_back(observed_snodes.size());

		if (monitor_type == TRACE) {
			m = new DensityTrace(observed_snodes, dim, density_type, monitor_name);
		}
		else if (monitor_type == MEAN) {
			m = new DensityMean(observed_snodes, dim, density_type, monitor_name);
		}
		else if (monitor_type == VARIANCE) {
			m = new DensityVariance(observed_snodes, dim, density_type, monitor_name);
		}
		else if (monitor_type == TOTAL) {
			m = new DensityTotal(observed_snodes, dim, density_type, monitor_name);
		}
		else if (monitor_type == POOLMEAN) {
			m = new DensityPoolMean(observed_snodes, dim, density_type, monitor_name);
		}
		else if (monitor_type == POOLVARIANCE) {
			m = new DensityPoolVariance(observed_snodes, dim, density_type, monitor_name);
		}
		else {
			throw std::logic_error("Unimplemented MonitorType in ObsStochDensMonitorFactory");
		}
		
		/* Set name attributes */
		
		// Will either be _observed_stochastic_ or deviance
		m->setName(name);

		// TOTAL is the only one summarised between variables
		if (monitor_type == TOTAL) {
		    m->setElementNames(vector<string>(1, monitor_name));
		}
		else {
			vector<string> const &onames = model->observedStochasticNodeNames();
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
