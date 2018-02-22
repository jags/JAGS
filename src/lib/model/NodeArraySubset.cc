#include <config.h>
#include <graph/Node.h>
#include <model/NodeArraySubset.h>
#include <model/NodeArray.h>
#include <sarray/RangeIterator.h>
#include <util/nainf.h>

#include <set>
#include <stdexcept>
#include <cmath>

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
	    _node_pointers = array->_node_pointers;
	    _offsets = array->_offsets;
	}
	else {
	    //Check validity of target range
	    if (!array->_range.contains(range)) {
		throw runtime_error(string("Cannot get subset ") +
				    array->_name + print(range) +
				    ". Range out of bounds");
	    }

	    for (RangeIterator p(range); !p.atEnd(); p.nextLeft()) {
		unsigned int i = array->_range.leftOffset(p);
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
	for (unsigned int i = 0; i < _node_pointers.size(); ++i) {
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
    
    vector<double> NodeArraySubset::logDensity(unsigned int chain, PDFType pdf_type, DensityType density_type) const
    {
	vector<double> ans;
	Node const *node = 0;
	double value = 0.0;
	for (unsigned int i = 0; i < _node_pointers.size(); ++i) {
	    if (_node_pointers[i]) {
		if (node != _node_pointers[i]) {
		    node = _node_pointers[i];
		    value = node->logDensity(chain, pdf_type);
			if ( density_type == DENSITY ) {
				value = std::exp(value);
			}
			else if ( density_type == DEVIANCE ) {
				value = -2.0 * value;
			}			
		}
		ans.push_back(value);
	    }
	    else {
		ans.push_back(JAGS_NA);
	    }
	}
	return ans;
    }
    
    vector<unsigned int> const &NodeArraySubset::dim() const
    {
	return _dim;
    }
    
    vector<Node const *> NodeArraySubset::nodes() const
    {
	vector<Node const *> ans;
	set<Node const *> nodeset;
	for (unsigned int i = 0; i < _node_pointers.size(); ++i) {
	    Node const * node = _node_pointers[i];
	    if (node && nodeset.insert(node).second) {
		ans.push_back(node);
	    }
	}
	return ans;
    }

    unsigned int NodeArraySubset::nchain() const
    {
	return _nchain;
    }

    unsigned int NodeArraySubset::length() const
    {
	return _node_pointers.size();
    }
    
    
} /* namespace jags */

