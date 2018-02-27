#include "DensityEnums.h"
#include "NodeDensityMonitorFactory.h"
#include "DensityTrace.h"
#include "DensityMean.h"
#include "DensityVariance.h"
#include "DensityTotal.h"
#include "DensityPoolMean.h"
#include "DensityPoolVariance.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <model/NodeArraySubset.h>
#include <sarray/RangeIterator.h>

#include <set>

using std::set;
using std::string;
using std::vector;

namespace jags {
namespace dic {

    Monitor *NodeDensityMonitorFactory::getMonitor(string const &name, 
						Range const &range,
						BUGSModel *model,
						string const &type,
						string &msg)
    {
		
		if (name == "_observed_stochastic_" ) {
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
		
		NodeArray *array = model->symtab().getVariable(name);
		if (!array) {
		    msg = string("Variable ") + name + " not found";
		    return 0;
		}
		
		/* Create the correct subtype of monitor */

		Monitor *m = 0;
		
		NodeArraySubset nodearray = NodeArraySubset(array, range);	
		
		if (monitor_type == TRACE) {
			m = new DensityTrace(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else if (monitor_type == MEAN) {
			m = new DensityMean(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else if (monitor_type == VARIANCE) {
			m = new DensityVariance(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else if (monitor_type == TOTAL) {
			m = new DensityTotal(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else if (monitor_type == POOLMEAN) {
			m = new DensityPoolMean(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else if (monitor_type == POOLVARIANCE) {
			m = new DensityPoolVariance(nodearray.allnodes(), nodearray.dim(), density_type, monitor_name);
		}
		else {
			throw std::logic_error("Unimplemented MonitorType in NodeDensityMonitorFactory");
		}
		
		/* Set name attributes */

		m->setName(name + print(range));
		Range node_range = range;
		if (isNULL(range)) {
		    //Special syntactic rule: a null range corresponds to the whole
		    //array
		    node_range = array->range();
		}

		// TOTAL is the only one summarised between variables
		if (monitor_type == TOTAL) {
		    m->setElementNames(vector<string>(1,type));
		}
		else {
			vector<string> elt_names;
			if (node_range.length() > 1) {
			    for (RangeIterator i(node_range); !i.atEnd(); i.nextLeft()) {
				elt_names.push_back(name + print(i));
			    }
			}
			else {
			    elt_names.push_back(name + print(range));
			}
			m->setElementNames(elt_names);
		}
		
		return m;
		
    }

    string NodeDensityMonitorFactory::name() const
    {
		return "dic::NodeDensity";
    }
	
}
}
