#include <config.h>
#include <compiler/ConstantFactory.h>
#include <graph/ConstantNode.h>
#include <compiler/NodeFactory.h>
#include <util/integer.h>

#include <sstream>
#include <cmath>
#include <cfloat>

using std::map;
using std::ostringstream;
using std::fabs;
using std::vector;
using std::sqrt;

ConstantFactory::ConstantFactory(unsigned int nchain)
    : _nchain(nchain)
{
}

ConstantNode * ConstantFactory::getConstantNode(double value, Graph &graph)
{
    ConstantNode *cnode = 0;

    map<double, ConstantNode*, ltdouble>::const_iterator i 
	= _constmap.find(value);
    if (i == _constmap.end()) {
	bool is_discrete;
	int ivalue = checkInteger(value, is_discrete);
	if (is_discrete) {
	    cnode = new ConstantNode(ivalue, _nchain);
	}
	else {
	    cnode = new ConstantNode(value, _nchain);
	}

	_constmap[value] = cnode;
	graph.add(cnode);
	return cnode;
    }
    else {
	return i->second;
    }
}

ConstantNode * ConstantFactory::getConstantNode(vector<unsigned int> const &dim,
						vector<double> const &value,
                                                Graph &graph)
{
    ConstantNode *cnode = 0;

    constpair cp(dim, value);

    map<constpair, ConstantNode*>::const_iterator i = _mv_constmap.find(cp);

    vector<double> ivalue = value;

    if (i == _mv_constmap.end()) {
	// Create a new constant node
	bool is_discrete = true;

	for (unsigned int j = 0; j < value.size(); ++j) {
	    ivalue[j] = checkInteger(value[j], is_discrete);
	    if (!is_discrete) {
		break;
	    }
	}
	
	if (is_discrete) {
	    cnode = new ConstantNode(dim, ivalue, _nchain);
	}
	else {
	    cnode = new ConstantNode(dim, value, _nchain);
	}
	
	_mv_constmap[cp] = cnode;
	graph.add(cnode);
	return cnode;
    }
    else {
	return i->second;
    }

}
