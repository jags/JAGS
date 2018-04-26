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

namespace jags {

    NodeArray::NodeArray(string const &name, vector<unsigned long> const &dim, 
			 unsigned int nchain)
	: _name(name), _range(dim), _nchain(nchain), 
	  _node_pointers(product(dim), 0),
	  _offsets(product(dim), numeric_limits<unsigned long>::max())
	  
    {
    }

    void NodeArray::insert(Node *node, Range const &target_range)
    {
	if (!node) {
	    throw logic_error(string("Attempt to insert NULL node at ") + 
			      name() + printRange(target_range));
	}
	if (node->dim() != target_range.dim(true)) {
	    throw runtime_error(string("Cannot insert node into ") + name() + 
				printRange(target_range) +
				". Dimension mismatch");
	}
	if (!_range.contains(target_range)) {
	    throw runtime_error(string("Cannot insert node into ") + name() + 
				printRange(target_range) +
				". Range out of bounds");
	}
	if (hasRepeats(target_range)) {
	    throw runtime_error(string("Cannot insert node into ") + name() +
				printRange(target_range) + 
				". Range has repeat indices");
	}

	/* Check that the range is not already occupied, even partially */
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft()) {
	    if (_node_pointers[_range.leftOffset(p)] != 0) {
		throw runtime_error(string("Node ") + name() 
				    + printRange(target_range)
				    + " overlaps previously defined nodes");
	    }
	}
	
	/* Set the _node_pointers array and the offset array */
	unsigned long k = 0;
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft())
	{
	    unsigned long i = _range.leftOffset(p);
	    _node_pointers[i] = node;
	    _offsets[i] = k++;
	}
	
	/* Add multivariate nodes to range map */
	if (node->length() > 1) {
	    _mv_nodes[target_range] = node;
	}
	
	/* Add node to the graph */
	_member_graph.insert(node);
    }
    
    Node *NodeArray::getSubset(Range const &target_range, Model &model)
    {
	//Check validity of target range
	if (!_range.contains(target_range)) {
	    throw runtime_error(string("Cannot get subset ") + name() + 
				printRange(target_range) +
				". Range out of bounds");
	}
	
	if (target_range.length() == 1) {
	    unsigned long start = _range.leftOffset(target_range.first());
	    Node *node = _node_pointers[start];
	    if (node && node->length() == 1) {
		if (_offsets[start] != 0) {
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
	    unsigned long i = _range.leftOffset(q);
	    if (_node_pointers[i] == 0) {
		return 0;
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
	if (!(_range == value.range())) {
	    throw runtime_error(string("Dimension mismatch in ") + name());
	}
	
	vector<double> const &x = value.value();
	unsigned long N = value.length();
	
	//Gather all the nodes for which a data value is supplied
	set<Node*> setnodes; 
	for (unsigned int i = 0; i < _range.length(); ++i) {
	    if (x[i] != JAGS_NA) {
		Node *node = _node_pointers[i];
		if (node == 0) {
		    string msg = "Attempt to set value of undefined node ";
		    throw runtime_error(msg + name() + 
					printIndex(value.range().leftIndex(i)));
		}
		switch(node->randomVariableStatus()) {
		case RV_FALSE:
		    throw NodeError(node, 
				    "Cannot set value of non-variable node");
		case RV_TRUE_OBSERVED:
		    throw NodeError(node,
				    "Cannot overwrite value of observed node");
		case RV_TRUE_UNOBSERVED:
		    setnodes.insert(node);		
		    break;
		}
	    }
	}
  

	for (set<Node*>::const_iterator p = setnodes.begin(); 
	     p != setnodes.end(); ++p) 
	{
	    //Step through each node
	    Node *node = *p;

	    vector<double> node_value(node->length());

	    //Get vector of values for this node
	    for (unsigned int i = 0; i < N; ++i) {
		if (_node_pointers[i] == node) {
		    if (_offsets[i] > node->length()) {
			throw logic_error("Invalid offset in NodeArray::setValue");
		    }
		    else {
			node_value[_offsets[i]] = x[i];
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

void NodeArray::getValue(SArray &value, unsigned int chain, 
			 bool (*condition)(Node const *)) const
{
    if (!(_range == value.range())) {
	string msg("Dimension mismatch when getting value of node array ");
	msg.append(name());
	throw runtime_error(msg);
    }

    unsigned long array_length = _range.length();
    vector<double> array_value(array_length);
    for (unsigned int j = 0; j < array_length; ++j) {
	Node const *node = _node_pointers[j];
	if (node && condition(node)) {
	    array_value[j] = node->value(chain)[_offsets[j]];
	}
	else {
	    array_value[j] = JAGS_NA;
	}
    }

    value.setValue(array_value);
}

void NodeArray::setData(SArray const &value, Model *model)
{
    if (!(_range == value.range())) {
	throw runtime_error(string("Dimension mismatch when setting value of node array ") + name());
    }

    vector<double> const &x = value.value();
  
    //Gather all the nodes for which a data value is supplied
    for (unsigned int i = 0; i < _range.length(); ++i) {
	if (x[i] != JAGS_NA) {
	    if (_node_pointers[i] == 0) {
		//Insert a new constant data node
		ConstantNode *cnode = new ConstantNode(x[i], _nchain, true);
		model->addNode(cnode);
		vector<unsigned long> index = _range.leftIndex(i);
		insert(cnode, SimpleRange(index, index));
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
	    for (unsigned int i = 0; i < _range.length(); ++i) {
		if (_node_pointers[i] == node) {
		    return SimpleRange(_range.leftIndex(i),
				       _range.leftIndex(i));
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
	
} //namespace jags
