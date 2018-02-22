#include "NodeDensityMonitorFactory.h"
#include "DensityTrace.h"
//#include "DensityMean.h"
//#include "DensityVariance.h"
//#include "DensityTotal.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
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
	
		/* Work out the precise type of monitor */
		
		// Used to help resolve aliases:
		string monitor_name = "";

		MonitorType monitor_type; 
		// enum declared in model/NodeArraySubset.h:
		DensityType density_type; 
		
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
		else {
			// If not listed above:
			return 0;
		}
		
		/*
		if (name != "deviance")
		    return 0;
		if (!isNULL(range)) {
		    msg = "cannot monitor a subset of deviance";
		    return 0;
	
		vector<StochasticNode *> const &snodes = model->stochasticNodes();
		vector<StochasticNode const *> observed_snodes;
		for (unsigned int i = 0; i < snodes.size(); ++i) {
		    if (snodes[i]->isFixed()) {
			observed_snodes.push_back(snodes[i]);
		    }
		}
		if (observed_snodes.empty()) {
		    msg = "There are no observed stochastic nodes";
		    return 0;
		}

		Monitor *m = 0;

		if (type == "mean") {
		    m = new DevianceMean(observed_snodes);
		    m->setName(name);
		    vector<string> onames(observed_snodes.size());
		    for (unsigned int i = 0; i < observed_snodes.size(); ++i) {
			onames[i] = model->symtab().getName(observed_snodes[i]);
		    }
		    m->setElementNames(onames);
		}
		else if (type == "trace") {
		    m = new DevianceTrace(observed_snodes);
		    m->setName("deviance");
		    m->setElementNames(vector<string>(1,"deviance"));
		}
		return m;

		}*/


		NodeArray *array = model->symtab().getVariable(name);
		if (!array) {
		    msg = string("Variable ") + name + " not found";
		    return 0;
		}

		/* Create the correct subtype of monitor */

		Monitor *m = 0;

		if (monitor_type == TRACE) {
			m = new DensityTrace(NodeArraySubset(array, range), density_type, monitor_name);
		}
/*		else if (monitor_type == MEAN) {
			m = new DensityMean(NodeArraySubset(array, range), density_type, monitor_name);
		}
		else if (monitor_type == VARIANCE) {
			m = new DensityVariance(NodeArraySubset(array, range), density_type, monitor_name);
		}
		else if (monitor_type == TOTAL) {
			m = new DensityTotal(NodeArraySubset(array, range), density_type, monitor_name);
		}
*/		else {
			throw std::logic_error("Unimplemented monitor_type in NodeDensityMonitorFactory");
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
	
}}
