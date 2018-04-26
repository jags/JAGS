#include <config.h>
#include <graph/Node.h>
#include <graph/NodeError.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <stdexcept>
#include <algorithm>

using std::string;
using std::vector;
using std::logic_error;
using std::copy;
using std::find;
using std::list;

namespace jags {

class DeterminsticNode;
class StochasticNode;

Node::Node(vector<unsigned long> const &dim, unsigned int nchain)
    : _parents(0), _stoch_children(0), _dtrm_children(0), 
      _dim(getUnique(dim)), _length(product(dim)), _nchain(nchain), _data(0)
{
    if (nchain==0)
	throw logic_error("Node must have at least one chain");

    unsigned long N = _length * _nchain;
    _data = new double[N];
    for (unsigned long i = 0; i < N; ++i) {
	_data[i] = JAGS_NA;
    }

    _dtrm_children = new list<DeterministicNode*>;
    _stoch_children = new list<StochasticNode*>;
}

Node::Node(vector<unsigned long> const &dim, unsigned int nchain,
	   vector<Node const *> const &parents)
    : _parents(parents), _stoch_children(0), _dtrm_children(0), 
      _dim(getUnique(dim)), _length(product(dim)),
      _nchain(nchain), _data(0)
{
    if (nchain==0)
	throw logic_error("Node must have at least one chain");

    unsigned long N = _length * _nchain;
    _data = new double[N];
    for (unsigned long i = 0; i < N; ++i) {
	_data[i] = JAGS_NA;
    }
  
    _stoch_children = new list<StochasticNode*>;
    _dtrm_children = new list<DeterministicNode*>;
}

Node::~Node()
{
    delete [] _data;
    delete _stoch_children;
    delete _dtrm_children;
}

vector <Node const *> const &Node::parents() const
{
    return _parents;
}

list<StochasticNode*> const *Node::stochasticChildren() 
{
    return _stoch_children;
}

list<DeterministicNode*> const *Node::deterministicChildren() 
{
    return _dtrm_children;
}

static bool isInitialized(Node const *node, unsigned int n)
{
    double const *value = node->value(n);
    for (unsigned long i = 0; i < node->length(); ++i) {
	if (value[i] == JAGS_NA) 
	    return false;
    }
    return true;
}

bool Node::initialize(RNG *rng, unsigned int n)
{
    // Test whether node is already initialized and, if so, skip it
    if (isInitialized(this, n))
        return true;

    // Check that parents are initialized
    for (unsigned long i = 0; i < _parents.size(); ++i) {
        if (!isInitialized(_parents[i], n)) {
	    return false; // Uninitialized parent
        }
    }

    randomSample(rng, n);

    return true; 
}
    
unsigned int Node::nchain() const
{
  return _nchain;
}

unsigned int countChains(vector<Node const *> const &parameters)
{
    unsigned int nchain = parameters.empty() ? 0 : parameters[0]->nchain();

    for (unsigned long i = 1; i < parameters.size(); ++i) {
	if (parameters[i]->nchain() != nchain) {
	    nchain = 0;
	    break;
	}
    }

    return nchain;
}

void Node::setValue(double const *value, unsigned long length, unsigned int chain)
{
   if (length != _length)
      throw NodeError(this, "Length mismatch in Node::setValue");
   if (chain >= _nchain)
      throw NodeError(this, "Invalid chain in Node::setValue");

   copy(value, value + _length, _data + chain * _length);
}

void Node::swapValue(unsigned int chain1, unsigned int chain2)
{
    double *value1 = _data + chain1 * _length;
    double *value2 = _data + chain2 * _length;
    for (unsigned int i = 0; i < _length; ++i) {
	double v = value1[i];
	value1[i] = value2[i];
	value2[i] = v;
    }
}

double const *Node::value(unsigned int chain) const
{
    return _data + chain * _length;
}

vector<unsigned long> const &Node::dim() const
{
    return _dim;
}

unsigned long Node::length() const
{
    return _length;
}

void Node::addChild(DeterministicNode *node) const
{
    _dtrm_children->push_back(node);
}

void Node::addChild(StochasticNode *node) const
{
    _stoch_children->push_back(node);
}

void Node::removeChild(DeterministicNode *node) const
{
    /* 
       Removes the given node from the list of deterministic children.
       
       When Model::~Model is called, all nodes are deleted in reverse
       order of construction. In this case, the element of _dtrm_node
       to remove is the last one. For efficiency, we therefore search
       for the element of _dtrm_node to remove starting at the
       end. (NB Searching from the beginning results in quadratic
       complexity in the size of _dtrm_node, which can cause real
       efficiency problems when deleting a model)
    */
    
    list<DeterministicNode*>::reverse_iterator p =
	find(_dtrm_children->rbegin(), _dtrm_children->rend(), node);
    if (p != _dtrm_children->rend()) {
	//The erase function only accepts iterators. Some shennanigans
	//are required to convert the reverse iterator correctly.
	_dtrm_children->erase((++p).base());
    }
}

void Node::removeChild(StochasticNode *node) const
{
    /* See comments in removeChild for DeterministicNodes */
    
    list<StochasticNode*>::reverse_iterator p = 
	find(_stoch_children->rbegin(), _stoch_children->rend(), node);
    if (p != _stoch_children->rend()) {
	_stoch_children->erase((++p).base());
    }
}

} //namespace jags
