#include <config.h>
#include <model/NodeArray.h>
#include <model/Model.h>
#include <graph/ConstantNode.h>
#include <graph/StochasticNode.h>
#include <graph/AggNode.h>
#include <sarray/RangeIterator.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <string>
#include <stdexcept>
#include <limits>

using std::pair;
using std::vector;
using std::map;
using std::string;
using std::runtime_error;
using std::logic_error;
using std::set;
using std::numeric_limits;

static bool hasRepeats(jags::Range const &target_range) 
{
    /* Returns true if the target range has any repeated indices 

       We choose the vectorized version of set::insert as it is
       amortized linear time in the length of the index vector
       scope[i] if the indices are in increasing order, which should
       be true most of the time.
    */
    
    vector<vector<unsigned long> > const &scope = target_range.scope();
    for (unsigned long i = 0; i < scope.size(); ++i) {
	set<unsigned long> seen;
	seen.insert(scope[i].begin(), scope[i].end());
	if (seen.size() != scope[i].size()) return true;
    }
    return false;
}

static vector<unsigned long> expand(vector<unsigned long> const &dim)
{
    vector<unsigned long> truedim = dim;
    for (unsigned long i = 0; i < dim.size(); ++i) {
	truedim[i] = static_cast<unsigned long>(1.05 * dim[i]);
    }

    return truedim;
}

namespace jags {
    
    NodeArray::NodeArray(string const &name, vector<unsigned long> const &dim, 
			 unsigned int nchain)
	: _name(name), _range(dim), _true_range(expand(dim)), _nchain(nchain), 
	  _node_pointers(product(expand(dim)), nullptr),
	  _offsets(product(expand(dim)), numeric_limits<unsigned long>::max()),
	  _locked(false)
	  
    {
    }

    void NodeArray::grow(Range const &target_range) {

	vector<unsigned long> upper = _range.upper();
	vector<unsigned long> true_upper = _true_range.upper();
	
	bool resize = false, extend = false;
	vector<vector<unsigned long> > const &scope = target_range.scope();
	for (unsigned long d = 0; d < scope.size(); ++d) {
	    vector<unsigned long> const &v = scope[d];
	    for (unsigned long i = 0; i < v.size(); ++i) {
		if (v[i] > true_upper[d]) {
		    upper[d] = v[i];
		    true_upper[d] = v[i];
		    resize = true;
		}
		else if (v[i] > upper[d]) {
		    upper[d] = v[i];
		    extend = true;
		}
	    }
	}

	if (resize) {
	    SimpleRange new_range = SimpleRange(upper);
	    SimpleRange new_true_range = SimpleRange(expand(upper));
	    unsigned long N  = product(new_true_range.dim(false));
	    vector<unsigned long> new_offsets(N, numeric_limits<unsigned long>::max());
	    vector<Node *> new_node_pointers(N, nullptr);
	    for (RangeIterator p(_range); !p.atEnd(); p.nextLeft()) {
		unsigned long k = _true_range.leftOffset(p);
		unsigned long l = new_true_range.leftOffset(p);
		new_offsets[l] = _offsets[k];
		new_node_pointers[l] = _node_pointers[k];
	    }

	    _range = new_range;
	    _true_range = new_true_range;
	    _node_pointers = new_node_pointers;
	    _offsets = new_offsets;
	}
	else if (extend) {
	    _range = SimpleRange(upper);
	}
	
    }
    
    void NodeArray::insert(Node *node, Range const &target_range)
    {
	// Check validity of target range
	if (hasRepeats(target_range)) {
	    throw runtime_error(string("Cannot insert node into ") + name() +
				printRange(target_range) + 
				". Range has repeat indices");
	}
	if (!_range.contains(target_range)) {
	    if (_locked) {
		throw runtime_error(string("Cannot insert node into ") + name()
				    + printRange(target_range) +
				    ". Range out of bounds");
	    }
	    else {
		grow(target_range);
	    }
	}
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft()) {
	    if (_node_pointers[_true_range.leftOffset(p)] != nullptr) {
		throw runtime_error(string("Node ") + name() 
				    + printRange(target_range)
				    + " overlaps previously defined nodes");
	    }
	}
	if (!node) {
	    //Was an error but now NodeArrays can grow dynamically it is
	    //useful to be able to pass a null node.
	    return;
	}
	
	if (node->dim() != target_range.dim(true)) {
	    throw runtime_error(string("Cannot insert node into ") + name() + 
				printRange(target_range) +
				". Dimension mismatch");
	}
	
	// Set the _node_pointers array and the offset array
	unsigned long s = 0;
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft())
	{
	    unsigned long k = _true_range.leftOffset(p);
	    _node_pointers[k] = node;
	    _offsets[k] = s++;
	}
	
	// Add multivariate nodes to range map
	if (node->length() > 1) {
	    _mv_nodes[target_range] = node;
	}
	
	// Add node to the graph
	_member_graph.insert(node);
    }
    
    Node *NodeArray::getSubset(Range const &target_range, Model &model)
    {
	//Check validity of target range
	if (!isNULL(target_range) && !_range.contains(target_range)) {
	    if (_locked) {
		throw runtime_error(string("Cannot get subset ") + name() + 
				    printRange(target_range) +
				    ". Range out of bounds");
	    }
	    else {
		return nullptr;
	    }
	}
	
	if (target_range.length() == 1) {
	    unsigned long i = _true_range.leftOffset(target_range.first());
	    Node *node = _node_pointers[i];
	    if (node && node->length() == 1) {
		if (_offsets[i] != 0) {
		    throw logic_error("Invalid scalar node in NodeArray");
		}
		return node;
	    }
	}
	else {
	    map<Range, Node *>::const_iterator p = _mv_nodes.find(target_range);
	    if (p != _mv_nodes.end()) {
		return p->second;
	    }
	}
	
	/* If range corresponds to a previously created subset, then
	 * return this */
	map<Range, AggNode *>::iterator p = _generated_nodes.find(target_range);
	if (p != _generated_nodes.end()) {
	    return p->second;
	}
	
	/* Otherwise create an aggregate node */
	
	vector<Node const *> nodes;
	vector<unsigned long> offsets;
	for (RangeIterator q(target_range); !q.atEnd(); q.nextLeft()) {
	    unsigned long i = _true_range.leftOffset(q);
	    if (_node_pointers[i] == nullptr) {
		return nullptr;
	    }
	    nodes.push_back(_node_pointers[i]);
	    offsets.push_back(_offsets[i]);
	}
	AggNode *anode = new AggNode(target_range.dim(true), _nchain,
				     nodes, offsets);
	_generated_nodes[target_range] = anode;
	model.addNode(anode);
	_member_graph.insert(anode);
	return anode;
    }

    void NodeArray::setValue(SArray const &value, unsigned int chain)
    {
	if (_range != value.range()) {
	    throw runtime_error(string("Dimension mismatch in ") + name());
	}
	
	vector<double> const &x = value.value();
	
	//Gather all the nodes for which a data value is supplied
	set<Node*> setnodes;
	for (RangeIterator p(_range); !p.atEnd(); p.nextLeft()) {
	    unsigned long j = _range.leftOffset(p);
	    if (x[j] != JAGS_NA) {
		unsigned long k = _true_range.leftOffset(p);
		Node *node = _node_pointers[k];
		if (node == nullptr) {
		    string msg = "Attempt to set value of undefined node ";
		    throw runtime_error(msg + name() + printIndex(p));
		}
		if (!node->isRandomVariable()) {
		    throw NodeError(node, 
				    "Cannot set value of non-variable node");
		}
		if (node->isObserved(_offsets[k])) {
		    throw NodeError(node,
				    "Cannot overwrite value of observed node");
		}
		setnodes.insert(node);		
	    }
	}
  

	for (set<Node*>::const_iterator p = setnodes.begin();
	     p != setnodes.end(); ++p) 
	{
	    //Step through each node
	    Node *node = *p;

	    vector<double> node_value(node->length());

	    //Get vector of values for this node
	    for (RangeIterator q(_range); !q.atEnd(); q.nextLeft()) {
		unsigned long k = _true_range.leftOffset(q);
		if (_node_pointers[k] == node) {
		    if (_offsets[k] > node->length()) {
			throw logic_error("Invalid offset in NodeArray::setValue");
		    }
		    else {
			unsigned long j = _range.leftOffset(q);
			node_value[_offsets[k]] = x[j];
		    }
		}
	    }
	    // If there are any missing values, they must all be missing
	    bool missing = node_value[0] == JAGS_NA;
	    for (unsigned int j = 1; j < node->length(); ++j) {
		if ((node_value[j] == JAGS_NA) != missing) {
		    throw NodeError(node,"Values supplied for node are partially missing");
		}
	    }
	    if (!missing) {
		node->setValue(&node_value[0], node->length(), chain);
	    }
	}
    }

    void
    NodeArray::getValue(SArray &value, unsigned int chain, ValueType type) const
    {
	if (_range != value.range()) {
	    string msg("Dimension mismatch when getting value of node array ");
	    msg.append(name());
	    throw runtime_error(msg);
	}

	vector<double> array_value(_range.length());
	for (RangeIterator p(_range); !p.atEnd(); p.nextLeft()) {
	    unsigned long j = _range.leftOffset(p);
	    unsigned long k = _true_range.leftOffset(p);
	    Node const *node = _node_pointers[k];
	    unsigned long offset = _offsets[k];
	    bool condition = false;
	    if (node != nullptr) {
		switch(type) {
		case DATA_VALUES:
		    condition = node->isRandomVariable() && node->isObserved(offset);
		break;
		case PARAMETER_VALUES:
		    condition = node->isRandomVariable() && !node->isObserved(offset);
		    break;
		case ALL_VALUES:
		    condition = true;
		    break;
		}
	    }
	    array_value[j] = condition ? node->value(chain)[offset] : JAGS_NA;
	}

	value.setValue(array_value);
    }

    void NodeArray::setData(SArray const &value, Model *model)
    {
	if (_range != value.range()) {
	    throw runtime_error(string("Dimension mismatch when setting value of node array ") + name());
	}

	vector<double> const &x = value.value();
  
	//Gather all the nodes for which a data value is supplied
	for (RangeIterator p(_range); !p.atEnd(); p.nextLeft()) {
	    unsigned long j = _range.leftOffset(p);
	    if (x[j] != JAGS_NA) {
		unsigned long k = _true_range.leftOffset(p);
		if (_node_pointers[k] == nullptr) {
		    //Insert a new constant data node
		    ConstantNode *cnode = new ConstantNode(x[j], _nchain, true);
		    model->addNode(cnode);
		    insert(cnode, SimpleRange(p, p));
		}
		else {
		    throw logic_error("Error in NodeArray::setData");
		}
	    }
	}
    }


    string const &NodeArray::name() const
    {
	return _name;
    }

    SimpleRange const &NodeArray::range() const
    {
	return _range;
    }

    Range NodeArray::getRange(Node const *node) const
    {
	if (!_member_graph.contains(node)) {
	    return Range();
	}
	
	//Look among inserted nodes first
	if (node->length() == 1) {
	    for (unsigned int i = 0; i < _true_range.length(); ++i) {
		if (_node_pointers[i] == node) {
		    return SimpleRange(_true_range.leftIndex(i),
				       _true_range.leftIndex(i));
		}
	    }
	}
	else {
	    for (map<Range, Node *>::const_iterator p = _mv_nodes.begin();
		 p != _mv_nodes.end(); ++p) 
	    { 
		if (node == p->second) {
		    return p->first;
		}
	    }
	}

	//Then among generated nodes
	map<Range, AggNode *>::const_iterator p = _generated_nodes.begin();
	for ( ; p != _generated_nodes.end(); ++p) 
	{ 
	    if (node == p->second) {
		break;
	    }
	}

	if (p == _generated_nodes.end()) {
	    throw logic_error("Failed to find Node range");
	}

	return p->first;
    }

    unsigned int NodeArray::nchain() const
    {
	return _nchain;
    }

    void NodeArray::lock()
    {
	_locked = true;
    }

    bool NodeArray::isLocked() const
    {
	return _locked;
    }
	
} //namespace jags
