#include <config.h>
#include <graph/Node.h>
#include <graph/NodeError.h>
#include <graph/NodeNameTab.h>
#include <sarray/SArray.h>

#include <stdexcept>

using std::set;
using std::string;
using std::vector;
using std::logic_error;

vector<SArray*> mkData(Index const &dim, unsigned int nchain)
{
  vector<SArray*> data(nchain);
  for (unsigned int n = 0; n < nchain; ++n) {
    data[n] = new SArray(dim);
  }
  return data;
}

Node::Node(Index const &dim, unsigned int nchain)
  : _parents(0), _children(), _ref(0), _data(mkData(dim,nchain))
{
  if (nchain==0)
    throw logic_error("Node must have at least one chain");
}

Node::Node(Index const &dim, vector<Node *> const &parents)
  : _parents(parents), _children(), _ref(0),
    _data(vector<SArray*>(mkData(dim,countChains(parents))))
{
  if (nchain() == 0) {
    throw logic_error("chain number mismatch in Node constructor");
  }

  for (unsigned int i = 0; i < parents.size(); ++i) {
    if (parents[i] == this) {
      throw NodeError(this, "Node cannot be its own parent");
    }
  }
  
  for (unsigned int i = 0; i < parents.size(); ++i) {
    parents[i]->_children.insert(this);
  }
}

Node::~Node()
{
  for (unsigned int n = 0; n < _data.size(); ++n) {
    delete _data[n];
  }

  for (unsigned int i = 0; i < _parents.size(); ++i) {
    _parents[i]->_children.erase(this);
  }
  for (set<Node*>::iterator p = _children.begin(); p != _children.end(); p++)
    {
      vector<Node*> &P = (*p)->_parents;
      vector<Node*>::iterator cp;
      while ( (cp = find(P.begin(), P.end(), this)) != P.end()) {
         P.erase(cp);
      }
    }
}

void Node::ref()
{
  _ref++;
}

void Node::unref()
{
  _ref--;
  if (_ref == 0) {
    delete this;
  }
}

unsigned int Node::refCount() const
{
  return _ref;
}

vector <Node *> const &Node::parents() const
{
  return _parents;
}


set<Node*> const &Node::children() const
{
  return _children;
}

static bool isInitialized(Node const *node, unsigned int chain)
{
   double const *value = node->value(chain);
   for (unsigned int i = 0; i < node->length(); ++i) {
      if (value[i] == JAGS_NA) {
          return false;
      }
   }
   return true;
}

bool Node::initialize()
{
    // Test whether node is already initialized and, if so, skip it
    bool isinitialized = true;
    for (unsigned int n = 0; n < nchain(); ++n) {
        if (!isInitialized(this, n)) {
            isinitialized = false;
            break;
        }
    }
    if (isinitialized) {
        return true; //Nothing to do
    }

    // Initialize each chain, if possible
    isinitialized = true;
    for (unsigned int n = 0; n < nchain(); ++n) {
        if (!isInitialized(this, n)) {
	    bool caninit = true;
	    for (unsigned int i = 0; i < _parents.size(); ++i) {
	        if (!isInitialized(_parents[i], n)) {
		    caninit = false;
		    break;
	        }
	    }
	    if (caninit) {
	        deterministicSample(n);
	    }
	    else {
	        isinitialized = false;
	    }
        }
    }

    if (!isinitialized)
	return false; //At least one chain uninitialized

    if (isVariable())
        return true; // Nothing more to do for variable nodes

    //Check if non-variable node is a function
    //of observed parents. In this case, it is observed also 
    for (unsigned int i = 0; i < _parents.size(); ++i) {
	if (!_parents[i]->isObserved()) {
            return true; //Not observed
	}
    }
    for (unsigned int n = 0; n < nchain(); ++n) {
        _data[n]->setFixed(true); //Observed
    }
    return true;
}
    
string Node::name(NodeNameTab const &name_table) const
{
  return name_table.getName(this);
}

/*
SArray const & Node::data(unsigned int chain) const
{
  return _data[chain];
}

SArray & Node::data(unsigned int chain)
{
  return _data[chain];
}
*/

Index const &Node::dim(bool drop) const
{
    return _data[0]->dim(drop);
}

unsigned long Node::length() const
{
    return _data[0]->length();
}

void Node::setObserved(double const *value, unsigned long length)
{
   for (unsigned int n = 0; n < nchain(); ++n) {
      _data[n]->setValue(value, length);
      _data[n]->setFixed(true);
   }
   //_status = NODE_DATA;
}

bool Node::isObserved() const
{
    /* A node is observed if all its data are fixed */
    for (unsigned int n = 0; n < nchain(); ++n) {
	if (!_data[n]->isFixed()) {
	    return false;
	}
    }
    return true;
}

unsigned int Node::nchain() const
{
  return _data.size();
}

unsigned int countChains(std::vector<Node*> const &parameters)
{
  unsigned int nchain = parameters[0]->nchain();

  for (unsigned int i = 1; i < parameters.size(); ++i) {
    if (parameters[i]->nchain() != nchain) {
      nchain = 0;
      break;
    }
  }

  return nchain;
}

void Node::setValue(double const *value, unsigned long length, unsigned int chain)
{
   _data[chain]->setValue(value, length);
}

bool Node::isDiscreteValued() const
{
   return _data[0]->isDiscreteValued();
}

void Node::setDiscreteValued()
{
   for (unsigned int n = 0; n < nchain(); ++n) {
      _data[n]->setDiscreteValued(true);
   }
}

double const *Node::value(unsigned int chain) const
{
   return _data[chain]->value();
}

SArray const *Node::data(unsigned int chain) const
{
  return _data[chain];
}
