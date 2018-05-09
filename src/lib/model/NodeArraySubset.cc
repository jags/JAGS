#include <config.h>
#include <graph/Node.h>
#include <model/NodeArraySubset.h>
#include <model/NodeArray.h>
#include <sarray/RangeIterator.h>
#include <util/nainf.h>

#include <set>
#include <stdexcept>
#include <sstream>


using std::set;
using std::vector;
using std::runtime_error;
using std::string;

namespace jags {

    NodeArraySubset::NodeArraySubset(NodeArray const *array, Range const &range)
	: _dim(range.dim(true)), _nchain(array->nchain())
    {
	if (isNULL(range)) {
	    // Special syntax rule: a NULL range means the whole array
	    _dim = array->range().dim(false);
	    for (RangeIterator p(array->_range); !p.atEnd(); p.nextLeft()) {
		unsigned long i = array->_true_range.leftOffset(p);
		_node_pointers.push_back(array->_node_pointers[i]);
		_offsets.push_back(array->_offsets[i]);
	    }
	}
	else {
		// Check that the implied number of dimensions is correct:
		unsigned long arraydim = array->_range.scope().size();
		unsigned long reqdim = range.scope().size();
		if ( arraydim != reqdim ) {
			std::ostringstream msgstr;
			msgstr << "Cannot get subset " << array->_name << 
			    printRange(range) <<
			    ": implied number of dimensions (" <<
			    reqdim <<
			    ") does not match the node dimensions (" <<
			    arraydim << ")";
			throw runtime_error(msgstr.str());
		}
		
	    //Check validity of target range
	    if (!array->_range.contains(range)) {
		throw runtime_error(string("Cannot get subset ") +
				    array->_name + printRange(range) +
				    ". Range out of bounds");
	    }

	    for (RangeIterator p(range); !p.atEnd(); p.nextLeft()) {
		unsigned long i = array->_true_range.leftOffset(p);
		_node_pointers.push_back(array->_node_pointers[i]);
		_offsets.push_back(array->_offsets[i]);
	    }
	}
    }
    
    vector<double> NodeArraySubset::value(unsigned int chain) const
    {
	vector<double> ans;
	Node const *node = 0;
	double const *values = 0;
	for (unsigned long i = 0; i < _node_pointers.size(); ++i) {
	    if (_node_pointers[i]) {
		if (node != _node_pointers[i]) {
		    node = _node_pointers[i];
		    values = node->value(chain);
		}
		ans.push_back(values[_offsets[i]]);
	    }
	    else {
		ans.push_back(JAGS_NA);
	    }
	}
	return ans;
    }
    
    vector<unsigned long> const &NodeArraySubset::dim() const
    {
	return _dim;
    }
    
    vector<Node const *> NodeArraySubset::nodes() const
    {
	vector<Node const *> ans;
	set<Node const *> nodeset;
	for (unsigned long i = 0; i < _node_pointers.size(); ++i) {
	    Node const * node = _node_pointers[i];
	    if (node && nodeset.insert(node).second) {
		ans.push_back(node);
	    }
	}
	return ans;
    }

    vector<Node const *> NodeArraySubset::allnodes() const
    {
	vector<Node const *> ans;
	for (unsigned int i = 0; i < _node_pointers.size(); ++i) {
	    Node const * node = _node_pointers[i];
		ans.push_back(node);
	}
	return ans;
    }

    unsigned int NodeArraySubset::nchain() const
    {
	return _nchain;
    }

    unsigned long NodeArraySubset::length() const
    {
	return _node_pointers.size();
    }
    
    
} /* namespace jags */

