#include <config.h>
#include <model/NodeArray.h>
#include <graph/ConstantNode.h>
#include <graph/AggNode.h>
#include <sarray/RangeIterator.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>

#include <string>
#include <stdexcept>

using std::vector;
using std::map;
using std::string;
using std::runtime_error;
using std::logic_error;
using std::set;

NodeArray::NodeArray(string const &name, Index const &dim, unsigned int nchain)
  : _name(name), _range(Index(dim.size()), dim), _nchain(nchain)
{
  long length = _range.length();
  _node_pointers = new Node *[length];
  _offsets = new long[length];
  for (long i = 0; i < length; i++) {
    _node_pointers[i] = 0;
    _offsets[i] = -1;
  }
}

NodeArray::~NodeArray()
{
  delete [] _node_pointers;
  delete [] _offsets;
}

bool NodeArray::isEmpty(Range const &target_range) const
{
  if (!_range.contains(target_range))
    throw logic_error("Range error in NodeArray::isEmpty");
  
  for (RangeIterator i(target_range); !i.atEnd(); i.nextLeft()) {
    if (_node_pointers[_range.leftOffset(i)] != 0)
      return false;
  }
  return true;
}

void NodeArray::insert(Node *node, Range const &target_range)
{
  if (!node) {
    throw logic_error(string("Attempt to insert NULL node at ") + name() +
		      print(target_range));
  }
  if (node->dim(true) != target_range.dim(true)) {
    throw runtime_error(string("Cannot insert node into ") + name() + 
			print(target_range) + ". Dimension mismatch");
  }
  if (!_range.contains(target_range)) {
    throw runtime_error(string("Cannot insert node into ") + name() + 
			print(target_range) + ". Range out of bounds");
  }
  if (!isEmpty(target_range)) {
    throw runtime_error(string("Node ") + name() + print(target_range)
			+ " overlaps previously defined nodes");
  }

  /* Set the _node_pointers array and the offset array */
  RangeIterator j(target_range);
  for (unsigned long k = 0; !j.atEnd(); j.nextLeft(), ++k)
    {
      long offset = _range.leftOffset(j);
      _node_pointers[offset] = node;
      _offsets[offset] = k;
    }

  /* Add to the graph */
  _graph.add(node);
}

Node *NodeArray::find(Range const &target_range) const
{
  if (!_range.contains(target_range)) {
    return 0;
  }

  long offset = _range.leftOffset(target_range.lower());
  Node *node = _node_pointers[offset];
  if (!node)
    return 0;

  if (node->dim(true) != target_range.dim(true))
    return 0;

  RangeIterator j(target_range);
  for (long k = 0; !j.atEnd(); j.nextLeft(), ++k) {
    offset = _range.leftOffset(j);
    if (_node_pointers[offset] != node || _offsets[offset] != k)
      return 0;
  }

  return node;
}

Node *NodeArray::getSubset(Range const &target_range)
{
  /* If range corresponds to a set node, then return this */
  Node *node = find(target_range);
  if (node)
    return node;

  /* If range corresponds to a previously created subset, then return this */
  map<Range, Node *>::iterator p = _generated_nodes.find(target_range);
  if (p != _generated_nodes.end()) {
    return p->second;
  }

  /* Otherwise create an aggregate node */
  vector<Node *> nodes;
  vector<unsigned long> offsets;
  for (RangeIterator i(target_range); !i.atEnd(); i.nextLeft()) {
    long offset = _range.leftOffset(i);
    if (_node_pointers[offset] == 0) {
      return 0;
    }
    nodes.push_back(_node_pointers[offset]);
    offsets.push_back(_offsets[offset]);
  }
  node = new AggNode(target_range.dim(false), nodes, offsets);
  _generated_nodes[target_range] = node;
  _graph.add(node);
  return node;
}

void NodeArray::setValue(SArray const &value, unsigned int chain,
			 set<Node*> &setnodes)

{
    if (!(_range == value.range())) {
	throw runtime_error(string("Dimension mismatch when setting value of node array ") + name());
    }
  
    double const *x = value.value();
    unsigned long N = value.length();

    //Gather all the nodes for which a data value is supplied
    for (unsigned long i = 0; i < _range.length(); ++i) {
	if (x[i] != JAGS_NA) {
	    Node *node = _node_pointers[i];
	    if (node == 0) {
		string msg = "Attempt to set value of undefined node";
		throw runtime_error(msg + name() + 
				    print(value.range().leftIndex(i)));
	    }
	    if (node->isVariable()) {
		setnodes.insert(node);
	    }
	    else {
		throw NodeError(node, 
				"Attempt to set value of non-variable node");
	    }
	}
    }
  
    set<Node*>::const_iterator p;
    double *node_value = new double[N];
    for (p = setnodes.begin(); p != setnodes.end(); ++p) {
	//Step through each node
	Node *node = *p;

	//Get vector of values for this node
	for (unsigned long i = 0; i < N; ++i) {
	    if (_node_pointers[i] == node) {
		if (_offsets[i] < 0 || _offsets[i] > node->length()) {
		    throw logic_error("Invalid offset in NodeArray::setValue");
		}
		else {
		    node_value[_offsets[i]] = x[i];
		}
	    }
	}
	// If there are any missing values, they must all be missing
	bool missing = node_value[0] == JAGS_NA;
	for (unsigned long j = 1; j < node->length(); ++j) {
	    if ((node_value[j] == JAGS_NA) != missing) {
		delete [] node_value;
		throw NodeError(node,"Values supplied for node are partially missing");
	    }
	}
	if (!missing) {
	    node->setValue(node_value, node->length(), chain);
	}
    }
    delete [] node_value;
}

void NodeArray::getValue(SArray &value, unsigned int chain, bool (*condition)(Node const *)) const
{
  if (!(_range == value.range())) {
    string msg("Dimension mismatch when getting value of node array ");
    msg.append(name());
    throw runtime_error(msg);
  }

  long array_length = _range.length();
  double *array_value = new double[array_length];
  for (int j = 0; j < array_length; ++j) {
    Node const *node = _node_pointers[j];
    if (node && condition(node)) {
      array_value[j] = node->value(chain)[_offsets[j]];
    }
    else {
      array_value[j] = JAGS_NA;
    }
  }


  value.setValue(array_value, array_length);
  delete [] array_value;
}

//FIXME: A lot of code overlap with setValue here.

#include <iostream>
void NodeArray::setData(SArray const &value)
{
  if (!(_range == value.range())) {
    throw runtime_error(string("Dimension mismatch when setting value of node array ") + name());
  }

  unsigned long N = value.length();  
  double const *x = value.value();
  
  //Gather all the nodes for which a data value is supplied
  set<Node*> setnodes;
  for (unsigned long i = 0; i < _range.length(); ++i) {
    if (x[i] != JAGS_NA) {
      if (_node_pointers[i] == 0) {
	//Insert a new constant node
	ConstantNode *cnode = new ConstantNode(x[i], _nchain);
	insert(cnode, _range.leftIndex(i));
      }
      else {
	//Existing node for which we must set value
	setnodes.insert(_node_pointers[i]);
      }
    }
  }
  
  set<Node*>::const_iterator p;
  double *node_value = new double[N];
  for (p = setnodes.begin(); p != setnodes.end(); ++p) {
    //Step through each node
    Node *node = *p;

    //Get vector of values for this node
    for (unsigned long i = 0; i < N; ++i) {
      if (_node_pointers[i] == node) {
	if (_offsets[i] < 0 || _offsets[i] > node->length()) {
	  throw logic_error("Invalid offset in NodeArray::setValue");
	}
	else {
	  node_value[_offsets[i]] = x[i];
	}
      }
    }
    // If there are any missing values, they must all be missing
    bool missing = (node_value[0] == JAGS_NA);
    for (unsigned long j = 1; j < node->length(); ++j) {
      if ((node_value[j] == JAGS_NA) != missing) {
	delete [] node_value;
	throw NodeError(node,"Values supplied for node are partially missing");
      }
    }
    if (!missing) {
      node->setObserved(node_value, node->length());
    }
  }
  delete [] node_value;
}

/* Old version

void NodeArray::setData(SArray const &value)
{
  if (!(_range == value.range())) {
    throw runtime_error(string("Dimension mismatch when setting data for node array ") + name());
  }
  
  double const *x = value.value();

  long j = 0;
  for (RangeIterator i(_range); !i.atEnd(); i.nextLeft(), ++j) {
    double fvalue = x[j];
    if (fvalue != JAGS_NA) {
      Node *node = _node_pointers[j];
      if (node == 0) {
	node = new ConstantNode(fvalue, nchain());
	insert(node, i);
      }
      else {
	for (unsigned int n = 0; n < nchain(); ++n) {
          if (!node->isObserved()) {
	    node->setValue(fvalue, _offsets[j], n);
	  }
          else {
            throw NodeError(node, "Attempt to change value of data node");
          }
        }
      }
    }
  }
  
  //Gather all the nodes for which a data value is supplied
  set<Node*> setnodes;
  for (unsigned long i = 0; i < _range.length(); ++i) {
    if (x[i] != JAGS_NA) {
      setnodes.insert(_node_pointers[i]);
    }
  }

  // Check to see if Node's value is partially missing
  for (set<Node*>::const_iterator p = setnodes.begin(); 
       p != setnodes.end(); ++p) 
    {
      double const *member_value = (*p)->data(0).value();
      unsigned long member_length = (*p)->data(0).length();
      for (unsigned long j = 0; j < member_length; ++j) {
	if (member_value[j] == JAGS_NA) {
	  throw NodeError(*p,"Values supplied for node are partially missing");
	}
      }     
      for (unsigned int n = 0; n < nchain(); ++n) {
	(*p)->data(n).setFixed(true);
      }
    }
}
*/

string const &NodeArray::name() const
{
  return _name;
}

Range const &NodeArray::range() const
{
  return _range;
}

Graph const &NodeArray::graph() const
{
  return _graph;
}

bool NodeArray::findActiveIndices(Index &ind, unsigned int k, 
				  Index const &lower, Index const &dim) const
{
  /* 
     We pay a heavy computational price for the flexibility of
     allowing users to insert multivariate nodes in arbritary
     ways into the NodeArray.
     
     Suppose we have an array of dimension [3,4,2,5], the lower index
     is [1,2,1,2] and the dimension of the node is [3,2]. Then the
     node could be inserted in 5 different ways.

     [1:3, 2:3, 1, 2] Active indices (0,1)
     [1:3, 2, 1:2, 2]                (0,2)
     [1:3, 2, 1, 2:3]                (0,3)
     [1, 2:4, 1:2, 2]                (1,2)
     [1, 2:4, 1, 2:3]                (1,3)

     We can't have active indices (2,3) because the node won't fit
  */

  if (k == 0)
    ind[k] = 0;
  else
    ind[k] = ind[k-1] + 1;
  unsigned int m = ind.size();
  unsigned int M = _range.ndim(false);
  for (;ind[k] + m <= M + k; ind[k] = ind[k] + 1) {
    if (k == m - 1) {
      Index upper(lower);
      for (unsigned int l = 0; l < m; ++l) {
	upper[ind[l]] = upper[ind[l]] + dim[l] - 1;
      }
      Range test_range(lower, upper);
      if (_range.contains(test_range)) {
	Node *node = _node_pointers[_range.leftOffset(lower)];
	int j = 0;
	bool ok = true;
	for (RangeIterator i(test_range); !i.atEnd(); i.nextLeft(), ++j) {
	  long offset = _range.leftOffset(i);
	  if (_node_pointers[offset] != node || _offsets[offset] != j) {
	    ok = false;
	    break;
	  }
	}
	if (ok)
	  return true;
      }
    }
    else {
      if (findActiveIndices(ind, k+1, lower, dim))
	return true;
    }
  }
  return false;
}

Range NodeArray::getRange(Node const *node) const
{
  if (!_graph.contains(node)) {
    return Range();
  }

  //Look in the generated nodes first
  for (map<Range, Node *>::const_iterator p = _generated_nodes.begin(); 
       p != _generated_nodes.end(); ++p) 
    { 
      if (node == p->second)
	return p->first;
    } 
  
  /* Find the lower limit of the range. This is easy */
  unsigned int ndim = _range.ndim(false);
  Index lower(ndim), upper(ndim);
  unsigned long j = 0;
  for (; j < _range.length(); ++j) {
    if (_node_pointers[j] == node) {
      lower = _range.leftIndex(j);
      break;
    }
  }
  if (j == _range.length()) {
    return Range();
  }

  unsigned int m = node->dim(false).size();
  Index ind(m);
  if (findActiveIndices(ind, 0, lower, node->dim(false))) {
    upper = lower;
    for (unsigned int l = 0; l < m; ++l) {
      upper[ind[l]] = upper[ind[l]] + node->dim(false)[l] - 1;
    }
    return Range(lower, upper);    
  }
  else {
    throw logic_error("Unable to find node range");
  }
}

unsigned int NodeArray::nchain() const
{
  return _nchain;
}
