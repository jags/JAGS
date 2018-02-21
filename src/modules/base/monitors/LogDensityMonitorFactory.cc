#include "LogDensityMonitorFactory.h"
#include "LogDensityMonitor.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <sarray/RangeIterator.h>

using std::set;
using std::string;
using std::vector;

namespace jags {
namespace base {

    Monitor *LogDensityMonitorFactory::getMonitor(string const &name,
					     Range const &range,
					     BUGSModel *model,
					     string const &type,
					     string &msg)
    {
	if (type != "logdensity")
	    return 0;

	NodeArray *array = model->symtab().getVariable(name);
	if (!array) {
	    msg = string("Variable ") + name + " not found";
	    return 0;
	}

	LogDensityMonitor *m = new LogDensityMonitor(NodeArraySubset(array, range));
	
	//Set name attributes 
	m->setName(name + print(range));
	Range node_range = range;
	if (isNULL(range)) {
	    //Special syntactic rule: a null range corresponds to the whole
	    //array
	    node_range = array->range();
	}
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
	
	return m;
    }

    string LogDensityMonitorFactory::name() const
    {
	return "base::LogDensity";
    }

}}
