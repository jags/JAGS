#include "DensityEnums.h"
#include "NodeDensityMonitorFactory.h"
#include "DensityTrace.h"
#include "DensityMean.h"
#include "DensityVariance.h"
#include "DensityTotal.h"
#include "DensityPoolMean.h"
#include "DensityPoolVariance.h"
#include "PenaltyPD.h"
#include "PenaltyPOPT.h"
#include "PenaltyPV.h"
#include "PenaltyPDTotal.h"
#include "PenaltyPOPTTotal.h"
#include "PenaltyPOPTTotalRep.h"

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
		
		// Should never be true but just in case:
		if (name == "_observations_" ) {
			return 0;
		}
		
		/* Work out the precise type of monitor */
				
		// enums declared in model/Monitor.h:
		MonitorType monitor_type = MTUNSET; 
		DensityType density_type = DTUNSET; 
		// Resolve the provided type:
		bool matched = getMonitorDensityTypes(type, monitor_type, density_type);
		if (!matched) {
			return 0;
		}
		
		/* Retrieve the node array  */
		
		NodeArray *array = model->symtab().getVariable(name);
		if (!array) {
		    msg = string("Variable ") + name + " not found";
		    return 0;
		}
		NodeArraySubset nodearray = NodeArraySubset(array, range);
		
		/* Do some checks and create the RNG vector for pd and popt monitors */

		vector<RNG*> rngs;
		if (monitor_type == PD || monitor_type == POPT
			|| monitor_type == PDTOTAL || monitor_type == POPTTOTAL ||
			monitor_type == POPTTOTALREP ) {

			if (model->nchain() < 2) {
			    msg = "at least two chains are required for a pD or popt monitor";
			    return 0;
			}
			
			/* 
			We could limit pD/popt monitors to observed stochastic nodes only
			But it does (maybe?) make sense as long as the parents of a node are unfixed
			Otherwise it comes out as 0 anyway - which makes sense (no parents are estimated)
			Note that pv can be calculated for any node with a density - which doesn't make sense
			if the parents are fixed
			TODO: create a node->areParentsFixed method to give an error for pv (and pD/popt??)
			
			// To limit pD / popt to observed stochastic nodes only (and pv if included above):
			vector<Node const *> const &reqnodes = nodearray.allnodes();
			for(unsigned int i = 0; i < reqnodes.size(); i++){
				if ( !reqnodes[i]->isStochastic() ) {
				    msg = "non-stochastic nodes cannot be included in an array subset for a pD or popt monitor";
				    return 0;
				}
				if ( !reqnodes[i]->isFixed() ) {
				    msg = "unobserved nodes cannot be included in an array subset for a pD or popt monitor";
				    return 0;
				}
			}*/
			
			for (unsigned int i = 0; i < model->nchain(); ++i) {
			    rngs.push_back(model->rng(i));
			}
		}

		/* Create the correct subtype of monitor */

		Monitor *m = 0;
		
		if (monitor_type == TRACE) {
			m = new DensityTrace(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == MEAN) {
			m = new DensityMean(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == VARIANCE) {
			m = new DensityVariance(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == TOTAL) {
			m = new DensityTotal(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == POOLMEAN) {
			m = new DensityPoolMean(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == POOLVARIANCE) {
			m = new DensityPoolVariance(nodearray.allnodes(), nodearray.dim(), density_type, type);
		}
		else if (monitor_type == PD) {
			m = new PenaltyPD(nodearray.allnodes(), nodearray.dim(), type, rngs, 10);
		}
		else if (monitor_type == POPT) {
			m = new PenaltyPOPT(nodearray.allnodes(), nodearray.dim(), type, rngs, 10);
		}
		else if (monitor_type == PDTOTAL) {
			m = new PenaltyPDTotal(nodearray.allnodes(), nodearray.dim(), type, rngs, 10);
		}
		else if (monitor_type == POPTTOTAL) {
			m = new PenaltyPOPTTotal(nodearray.allnodes(), nodearray.dim(), type, rngs, 10);
		}
		else if (monitor_type == POPTTOTALREP) {
			m = new PenaltyPOPTTotalRep(nodearray.allnodes(), nodearray.dim(), type, rngs, 10);
		}
		else if (monitor_type == PV) {
			m = new PenaltyPV(nodearray.allnodes(), nodearray.dim(), type);
		}
		else {
			throw std::logic_error("Unimplemented MonitorType in NodeDensityMonitorFactory");
		}
		
		if (!m ) {
			return m;
		}
		
		/* Set name attributes */

		m->setName(name + print(range));
		Range node_range = range;
		if (isNULL(range)) {
		    //Special syntactic rule: a null range corresponds to the whole
		    //array
		    node_range = array->range();
		}

		// These types are summarised between variables:
		if (monitor_type == TOTAL || monitor_type == PDTOTAL
			|| monitor_type == POPTTOTAL || monitor_type == PV) {
		    m->setElementNames(vector<string>(1, type));
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
