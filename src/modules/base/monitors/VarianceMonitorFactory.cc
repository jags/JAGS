#include "VarianceMonitorFactory.h"
#include "VarianceMonitor.h"
#include "PoolVarianceMonitor.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <sarray/RangeIterator.h>

using std::set;
using std::string;
using std::vector;

namespace jags {
namespace base {

    Monitor *VarianceMonitorFactory::getMonitor(string const &name,
					     Range const &range,
					     BUGSModel *model,
					     string const &type,
					     string &msg)
    {
	if (type != "variance" && type != "poolvariance")
	    return 0;

	NodeArray *array = model->symtab().getVariable(name);
	if (!array) {
	    msg = string("Variable ") + name + " not found";
	    return 0;
	}
	
	Monitor *m = 0;
	
	if ( type == "variance" ) {
		m = new VarianceMonitor(NodeArraySubset(array, range));
	}
	else if ( type == "poolvariance" ) {
		m = new PoolVarianceMonitor(NodeArraySubset(array, range));
	}
	else {
		throw std::logic_error("Unimplemented MonitorType in VarianceMonitorFactory");
	}
	
	//Set name attributes 
	m->setName(name + printRange(range));
	Range node_range = range;
	if (isNULL(range)) {
	    //Special syntactic rule: a null range corresponds to the whole
	    //array
	    node_range = array->range();
	}
	vector<string> elt_names;
	if (node_range.length() > 1) {
	    for (RangeIterator i(node_range); !i.atEnd(); i.nextLeft()) {
		elt_names.push_back(name + printIndex(i));
	    }
	}
	else {
	    elt_names.push_back(name + printRange(range));
	}
	m->setElementNames(elt_names);
	
	return m;
    }

    string VarianceMonitorFactory::name() const
    {
	return "base::Variance";
    }

}}
