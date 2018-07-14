#include "TraceMonitorFactory.h"
#include "TraceMonitor.h"

#include <model/BUGSModel.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <sarray/RangeIterator.h>

using std::set;
using std::string;
using std::vector;

namespace jags {
namespace base {

    Monitor *TraceMonitorFactory::getMonitor(string const &name,
					     Range const &range,
					     BUGSModel *model,
					     string const &type,
					     string &msg)
    {
	if (type != "trace")
	    return nullptr;

	NodeArray *array = model->symtab().getVariable(name);
	if (!array) {
	    msg = string("Variable ") + name + " not found";
	    return nullptr;
	}

	TraceMonitor *m = new TraceMonitor(NodeArraySubset(array, range));
	
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

    string TraceMonitorFactory::name() const
    {
	return "base::Trace";
    }

}}
