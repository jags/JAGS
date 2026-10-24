#include <config.h>
#include <model/NodeArray.h>
#include <model/Model.h>
#include <graph/ConstantNode.h>
#include <graph/StochasticNode.h>
#include <graph/AggNode.h>
#include <graph/MixtureNode.h>
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
    /* 
       Add some padding around the array dimensions to allow expansion
       without reallocating memory. This only affects larger arrays.
    */
    vector<unsigned long> truedim = dim;
    for (unsigned long i = 0; i < dim.size(); ++i) {
	truedim[i] += dim[i]/20;
    }

    return truedim;
}

namespace jags {
    
    /* 
       Given a vector of subset indices, this function modifies the argument
       "subsets", so that it contains a vector of pairs:
       - Index: Possible value of the variable indices
       - Range: Corresponding range of the array to take as a subset
       
       All possible values of Index are included in the vector, based on the
       the range of the array we are taking subsets from (default_range).
    */
    static void getSubsetRanges(vector<pair<vector<unsigned long>, Range>> &subsets,  
				vector<StochasticIndex> const &limits,
				SimpleRange const &default_range)
    {
	unsigned long ndim = limits.size();

	// Create upper and lower bounds
	vector<unsigned long> var_offset;
	vector<unsigned long> var_lower, var_upper;
	vector<vector<unsigned long> > scope(ndim);
	for (unsigned int j = 0; j < ndim; ++j) {
	    if (limits[j].isVariable()) {
		var_offset.push_back(j);
		var_lower.push_back(default_range.lower()[j]);
		var_upper.push_back(default_range.upper()[j]);
	    }
	    else {
		scope[j] = limits[j].fixedIndex();
	    }
	}

	SimpleRange var_range(var_lower, var_upper); //range of variable indices
	for (RangeIterator p(var_range); !p.atEnd(); p.nextLeft()) {
	    for (unsigned int k = 0; k < var_offset.size(); ++k) {
		scope[var_offset[k]] = vector<unsigned long>(1, p[k]);
	    }
	    subsets.push_back(pair<vector<unsigned long>, Range>(p, Range(scope)));
	}
    }

    NodeArray::NodeArray(string const &name, vector<unsigned long> const &dim, 
			 unsigned int nchain)
	: _name(name), _range(dim), _true_range(expand(dim)), _nchain(nchain), 
	  _node_pointers(_true_range.length(), nullptr),
	  _offsets(_true_range.length(), numeric_limits<unsigned long>::max()),
	  _locked(false)
    {
    }

    void NodeArray::grow(Range const &target_range) {

	if (_range.ndim(false) != target_range.ndim(false)) {
	    throw logic_error("Incorrect dimensions in NodeArray::grow");
	}
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
	    unsigned long N  = new_true_range.length();
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
	if (target_range.ndim(false) != _range.ndim(false)) {
	    throw runtime_error(string("Cannot insert node into ") + name() +
				printRange(target_range) +
				". Incorrect number of dimensions");
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
	if (!isNULL(target_range)) {
	    if (_range.ndim(false) != target_range.ndim(false)) {
		throw runtime_error(string("Cannot get subset ") + name() +
				    printRange(target_range) +
				    ". Incorrect number of dimensions");
	    }
	    if (!_range.contains(target_range)) {
		if (_locked) {
		    throw runtime_error(string("Cannot get subset ") + name() + 
					printRange(target_range) +
					". Range out of bounds");
		}
		else {
		    return nullptr;
		}
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

    Node * NodeArray::getMixture(vector<StochasticIndex> const &indices,
				 Model &model)
    {
	if (!_locked) {
	    //Can't create MixtureNodes until we stop growing the NodeArray
	    return nullptr;
	}
	
	//See if the mixture node exists already
	auto p = _mixture_nodes.find(indices);
	if (p != _mixture_nodes.end()) {
	    return p->second;
	}

	//Separate fixed and variable indices
	vector<vector<unsigned long>> fixed_indices, full_indices;
	vector<Node const *> variable_indices;
	for (unsigned int i = 0; i < indices.size(); ++i) {
	    if (!indices[i].isVariable()) {
		fixed_indices.push_back(indices[i].fixedIndex());
		full_indices.push_back(indices[i].fixedIndex());
	    }
	    else {
		variable_indices.push_back(indices[i].variableIndex());
		full_indices.push_back(vector<unsigned long>());
	    }
	}

	//All MixtureNodes with the same fixed indices can share the
	//same MixTab. See if a suitable MixTab exists already.
	auto q = _mix_tabs.find(full_indices);
	if (q == _mix_tabs.end()) {
	    //No pre-existing MixTab. So create one
	    vector<pair<vector<unsigned long>, Range>> subsets;  
	    getSubsetRanges(subsets, indices, this->range());

	    map<vector<unsigned long>, Node const *> mixmap;
	    Node *subset0 = nullptr; // Checking for degenerate mixture
	    for (unsigned int i = 0; i < subsets.size(); ++i) {
		Node *subset_node =
		    this->getSubset(subsets[i].second, model);
		if (subset_node) {
		    if (i == 0) {
			subset0 = subset_node;
		    }
		    else if (subset_node != subset0) {
			subset0 = nullptr;
		    }
		    mixmap[subsets[i].first] = subset_node;
		}
		else {
		    /* FIXME: In debug mode we'll have to find
		       something smarter
		    */
		    return nullptr;
		}
	    }
	    if (subset0 != nullptr) {
		// Degenerate mixture with all subsets equal to subset0.
		// No need to create a MixtureNode
		return subset0;
	    }
	    MixTab table(mixmap);
	    _mix_tabs.insert({full_indices, table});
	    q = _mix_tabs.find(full_indices);
	}
	
	//Create new MixtureNode and add it to the Model
	MixtureNode * mixnode = new MixtureNode(indices, _nchain, q->second);
	model.addNode(mixnode);
	
	//Insert into map for future calls to getMixture
	_mixture_nodes[indices] = mixnode;

	return mixnode;
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
	    if (!jags_isna(x[j])) {
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
				    "Cannot overwrite observed value");
		}
		setnodes.insert(node);		
	    }
	}

	for (set<Node*>::const_iterator p = setnodes.begin();
	     p != setnodes.end(); ++p) 
	{
	    //Step through each node
	    Node *node = *p;

	    /* Start with a copy of the current node value. We need to
	       do this when initializing partly observed nodes as we have
	       a mixture of already set data values and newly set initial
	       values */
	    vector<double> node_value(node->length());
	    double const *v = node->value(chain);
	    copy(v, v + node->length(), node_value.begin());
	    
	    //Get vector of values for this node
	    for (RangeIterator q(_range); !q.atEnd(); q.nextLeft()) {
		unsigned long k = _true_range.leftOffset(q);
		if (_node_pointers[k] == node) {
		    if (_offsets[k] > node->length()) {
			throw logic_error("Invalid offset in NodeArray::setValue");
		    }
		    else {
			//Only overwrite values if RHS is not missing
			unsigned long j = _range.leftOffset(q);
			if (!jags_isna(x[j])) {
			    node_value[_offsets[k]] = x[j];
			}
		    }
		}
	    }
	    node->setValue(&node_value[0], node->length(), chain);
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
	    if (!jags_isna(x[j])) {
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

    vector<StochasticIndex>
    NodeArray::getStochasticRange(Node const *node) const
    {
	for (auto p = _mixture_nodes.begin(); p != _mixture_nodes.end(); ++p) {
	    if (p->second == node)
		return p->first;
	}
	return vector<StochasticIndex>();
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
