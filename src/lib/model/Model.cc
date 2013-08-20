#include <config.h>
#include <model/Model.h>
#include <model/MonitorFactory.h>
#include <model/Monitor.h>
#include <sampler/Sampler.h>
#include <sampler/SamplerFactory.h>
#include <rng/RNGFactory.h>
#include <rng/RNG.h>
#include <graph/GraphMarks.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/ConstantNode.h>
#include <graph/NodeError.h>
#include <graph/Node.h>
#include <util/nainf.h>

#include <fstream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <functional>
#include <map>

using std::map;
using std::pair;
using std::binary_function;
using std::sort;
using std::vector;
using std::list;
using std::set;
using std::ofstream;
using std::logic_error;
using std::runtime_error;
using std::string;
using std::ostringstream;
using std::stable_sort;
using std::copy;
using std::min;
using std::max;
using std::reverse;

Model::Model(unsigned int nchain)
    : _samplers(0), _nchain(nchain), _rng(nchain, 0), _iteration(0),
      _is_initialized(false), _adapt(false), _data_gen(false)
{
}

Model::~Model()
{
    while(!_samplers.empty()) {
	Sampler *sampler0 = _samplers.back();
	delete sampler0;
	_samplers.pop_back();
    }
    for (list<Monitor*>::const_iterator p = _default_monitors.begin();
	 p != _default_monitors.end(); ++p)
    {
	delete *p;
    }

    //Delete nodes in reverse sampling order
    vector<Node*> managed_nodes;
    _graph.getSortedNodes(managed_nodes);
    while(!managed_nodes.empty())
    {
	Node *node = managed_nodes.back();
	delete node;
	managed_nodes.pop_back();
    }

}

Graph const &Model::graph() 
{
    return _graph;
}

bool Model::isInitialized()
{
    return _is_initialized;
}

void Model::chooseRNGs()
{
    /* Assign default RNG objects for any chain that does not
       currently have one */

    //Count number of unassigned RNGs
    unsigned int n  = 0;
    for (unsigned int i = 0; i < _nchain; ++i) {
	if (_rng[i] == 0)
	    ++n;
    }

    vector<RNG*> new_rngs;
    list<pair<RNGFactory*, bool> >::const_iterator p;
    for (p = rngFactories().begin(); p != rngFactories().end(); ++p) {
	if (p->second) {
	    vector<RNG*> rngs = p->first->makeRNGs(n);
	    if (rngs.size() > n) {
		throw logic_error("Too many rngs produced by RNG factory");
	    }
	    else {
		n -= rngs.size();
	    }
	    for (unsigned int j = 0; j < rngs.size(); ++j) {
		new_rngs.push_back(rngs[j]);
	    }
	    if (n == 0)
		break;
	}
    }
    
    if (n > 0) {
	throw runtime_error("Cannot generate sufficient RNGs");
    }
    else {
	unsigned int j = 0;
	for (unsigned int i = 0; i < _nchain; ++i) {
	    if (_rng[i] == 0) {
		_rng[i] = new_rngs[j++];
	    }
	}
    }
}

void Model::initialize(bool datagen)
{
    if (_is_initialized)
	throw logic_error("Model already initialized");

    if (!_graph.isClosed())
	throw runtime_error("Graph not closed");

    // Choose random number generators
    chooseRNGs();

    //Initialize nodes
    initializeNodes();
    
    // Check initial values of all stochastic nodes
    // Note that we need to do this before choosing samplers.
    if (!datagen) {
	for (unsigned int ch = 0; ch < _nchain; ++ch) {
	    for (unsigned int i = 0; i < _stochastic_nodes.size(); ++i) {
		StochasticNode const *snode = _stochastic_nodes[i];
		double ld = snode->logDensity(ch, PDF_PRIOR);
		if (jags_isnan(ld)) {
		    string msg = "Error calculating log density at initial values";
		    throw NodeError(snode, msg);
		}
		else if (ld == JAGS_NEGINF || (!jags_finite(ld) && ld < 0)) {
		    string msg;
		    if (snode->isObserved()) {
			msg = "Observed node";
		    }
		    else {
			msg = "Unobserved node";
		    }
		    msg.append(" inconsistent with ");
		    bool obs_par = true;
		    for (unsigned int j = 0; j < snode->parents().size(); ++j) {
			if (!snode->parents()[j]->isObserved()) {
			    obs_par = false;
			    break;
			}
		    }
		    if (obs_par) {
			msg.append("observed parents");
		    }
		    else {
			msg.append("unobserved parents");
		    }
		    msg.append(" at initialization.\n");
		    msg.append("Try setting appropriate initial values.");
		    throw NodeError(snode, msg);
		}
	    }
	}
    }

    // Choose Samplers
    chooseSamplers();
    
    if (datagen) {
	Graph egraph;
	for (set<Node *>::const_iterator p = _extra_nodes.begin(); 
	     p != _extra_nodes.end(); ++p)
	{
	    egraph.add(*p);
	}
	_sampled_extra.clear();
	egraph.getSortedNodes(_sampled_extra);
	_data_gen = true;
    }

    // Switch to adaptive mode if we find an adaptive sampler
    for (unsigned int i = 0; i < _samplers.size(); ++i) {
	if (_samplers[i]->isAdaptive()) {
	    _adapt = true;
	    break;
	}
    }
    
    _is_initialized = true;
}

void Model::initializeNodes()
{
    //Get nodes in forward-sampling order
    vector<Node*> sorted_nodes;
    _graph.getSortedNodes(sorted_nodes);

    vector<Node*>::const_iterator i;
    for (i = sorted_nodes.begin(); i != sorted_nodes.end(); ++i) {
	Node *node = *i;
	for (unsigned int n = 0; n < _nchain; ++n) {
	    if (!node->checkParentValues(n)) {
		throw NodeError(node, "Invalid parent values");
	    }
	    if (!node->initialize(n)) {
		throw NodeError(node, "Initialization failure");
	    } 
	}
    }
}

struct less_sampler {  
    /* 
       Comparison operator for Samplers which sorts them in 
       order according to the supplied sampler map
    */
    map<Sampler const*, unsigned int>  const & _sampler_map;

    less_sampler(map<Sampler const*, unsigned int> const &sampler_map) 
	: _sampler_map(sampler_map) {};

    bool operator()(Sampler const *x, Sampler const *y) const {
	return _sampler_map.find(x)->second < _sampler_map.find(y)->second;
    };

};

void Model::chooseSamplers()
{
    /*
     * Selects samplers. Samplers are selected by traversing the list
     * of SamplerFactories in order. If there are any informative
     * stochastic nodes left without samplers after all factories have
     * been tried, then a runtime error is thrown
     *
     * @see Model#samplerFactories
     */

    GraphMarks marks(_graph);
    Graph sample_graph;

    // Add observed stochastic nodes to the sample graph and mark
    // the informative nodes

    vector<Node const*> informative;
    vector<StochasticNode*>::const_iterator p;
    for (p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	if ((*p)->isObserved()) {
	    sample_graph.add(*p);
	    informative.push_back(*p);
	}
    }
    marks.markAncestors(informative, 1);

    for (p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	if ((*p)->isObserved()) {
	    marks.mark(*p, 2);
	}
    }

    //Triage on marked nodes. We do this twice: once for stochastic
    //nodes and once for all nodes.

    set<StochasticNode*> sset;
    for(p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	switch(marks.mark(*p)) {
	case 0:
	    _extra_nodes.insert(*p);
	    break;
	case 1:
	    sset.insert(*p); 
	    sample_graph.add(*p);
	    break;
	case 2:
	    sample_graph.add(*p);
	    break;
	default:
	    throw logic_error("Invalid mark");
	}
    }

    set<Node*> const &nodes = _graph.nodes();
    set<Node*>::const_iterator j;
    for (j = nodes.begin(); j != nodes.end(); ++j) {
	switch(marks.mark(*j)) {
	case 0:
	    _extra_nodes.insert(*j);
	    break;
	case 1: case 2:
	    sample_graph.add(*j);
	    break;
	}
    }

    // Traverse the list of samplers, selecting nodes that can be sampled
    list<pair<SamplerFactory *, bool> > const &sf = samplerFactories();
    for(list<pair<SamplerFactory *, bool> >::const_iterator q = sf.begin();
	q != sf.end(); ++q) 
    {
	if (q->second) {
	    vector<Sampler*> svec = q->first->makeSamplers(sset, sample_graph);
	    while (!svec.empty()) {
		for (unsigned int i = 0; i < svec.size(); ++i) {
		    vector<StochasticNode*> const &nodes = svec[i]->nodes();
		    for (unsigned int j = 0; j < nodes.size(); ++j) {
			sset.erase(nodes[j]);
		    }
		    _samplers.push_back(svec[i]);
		}
		svec = q->first->makeSamplers(sset, sample_graph);
	    }
	}
    }
  
    // Make sure we found a sampler for all the nodes
    if (!sset.empty()) {
	throw NodeError(*sset.begin(),
			"Unable to find appropriate sampler");
    }
  
    // Samplers are sorted in reverse sampling order: i.e. samplers
    // that are closer to the data are updated before samplers that
    // only affect higher-order parameters
    
    // Create a map associating each stochastic node with its index
    // in the vector _stochastic_nodes, corresponding to the order
    // in which they were added to the model
    map<StochasticNode const *, unsigned int> snode_map;
    for (unsigned int i = 0; i < _stochastic_nodes.size(); ++i) {
	snode_map[_stochastic_nodes[i]] = i;
    }

    // Create a map associating each sampler with the minimal index
    // of its sampled nodes.
    map<Sampler const *, unsigned int> sampler_map;
    for (unsigned int i = 0; i < _samplers.size(); ++i) {
	unsigned int min_index = _stochastic_nodes.size();
	vector<StochasticNode*> const &snodes = _samplers[i]->nodes();
	for (unsigned int j = 0; j < snodes.size(); ++j) {
	    map<StochasticNode const*, unsigned int>::const_iterator p 
		= snode_map.find(snodes[j]);
	    if (p == snode_map.end()) {
		throw logic_error("Invalid stochastic node map");
	    }
	    if (p->second < min_index) {
		min_index = p->second;
	    }
	}
	sampler_map[_samplers[i]] = min_index;
    }

    stable_sort(_samplers.begin(), _samplers.end(), less_sampler(sampler_map));
    reverse(_samplers.begin(), _samplers.end());
}

void Model::update(unsigned int niter)
{
    if (!_is_initialized) {
	throw logic_error("Attempt to update uninitialized model");
    }

    list<MonitorControl>::iterator p;
    for (p = _monitors.begin(); p != _monitors.end(); ++p) {
	p->reserve(niter);
    }

    for (unsigned int iter = 0; iter < niter; ++iter) {    
	
	for (vector<Sampler*>::iterator i = _samplers.begin(); 
	     i != _samplers.end(); ++i) 
	{
	    (*i)->update(_rng);
	}

	for (unsigned int n = 0; n < _nchain; ++n) {
	    for (vector<Node*>::const_iterator k = _sampled_extra.begin();
		 k != _sampled_extra.end(); ++k)
	    {
                if (!(*k)->checkParentValues(n)) {
                     throw NodeError(*k, "Invalid parent values");
                }
		(*k)->randomSample(_rng[n], n);
	    }
	}
	_iteration++;

	for (list<MonitorControl>::iterator k = _monitors.begin(); 
	     k != _monitors.end(); k++) 
	{
	    k->update(_iteration);
	}
    }
}

unsigned int Model::iteration() const
{
  return _iteration;
}

void Model::adaptOff() 
{
    for (vector<Sampler*>::const_iterator p = _samplers.begin(); 
	 p != _samplers.end(); ++p)
    {
	(*p)->adaptOff();
    }
    _adapt = false;
}

bool Model::checkAdaptation() const
{
    for (vector<Sampler*>::const_iterator p = _samplers.begin(); 
	 p != _samplers.end(); ++p)
    {
	if (!(*p)->checkAdaptation()) return false;
    }
    return true;
}

bool Model::isAdapting() const
{
  return _adapt;
}


void Model::setSampledExtra()
{
    /* If a mode is not a data generating model, uninformative nodes
       do not need to be updated, unless they have a descendant that 
       is being monitored. This function finds those nodes and adds
       them to the vector _sampled_extra.
    */
       
    if (_data_gen) {
	// In a data generating model, all uninformative nodes are
	// sampled, so nothing to be done
	return;
    }

    // Recalculate the vector of uninformative nodes that need sampling
	
    //Insert extra nodes into a new graph
    Graph egraph;
    for (set<Node *>::const_iterator p = _extra_nodes.begin(); 
	 p != _extra_nodes.end(); ++p)
    {
	egraph.add(*p);
    }
    //Mark the ancestors of all monitored nodes in this graph
    GraphMarks emarks(egraph);
    for (list<MonitorControl>::const_iterator p = _monitors.begin();
	 p != _monitors.end(); ++p)
    {
	for (unsigned int i = 0; i < p->monitor()->nodes().size(); ++i) {
	    Node const *node = p->monitor()->nodes()[i];
	    if (egraph.contains(node)) {
		emarks.mark(node, 1);
		//FIXME: call once
		emarks.markAncestors(vector<Node const *>(1, node), 1);
	    }
	}
    }
    //Remove unmarked nodes from graph
    for (set<Node *>::const_iterator p = _extra_nodes.begin(); 
	 p != _extra_nodes.end(); ++p)
    {
	if (emarks.mark(*p) == 0)
	    egraph.remove(*p);
    }
    //Replace vector of sampled extra nodes
    _sampled_extra.clear();
    egraph.getSortedNodes(_sampled_extra);
}

void Model::addMonitor(Monitor *monitor, unsigned int thin)
{
    if (_adapt) {
	throw runtime_error("Turn off adaptive mode before setting monitors");
    }
    
    _monitors.push_back(MonitorControl(monitor, _iteration+1, thin));
    setSampledExtra();
}

void Model::removeMonitor(Monitor *monitor)
{
    for(list<MonitorControl>::iterator p = _monitors.begin();
	p != _monitors.end(); ++p)
    {
	if (p->monitor() == monitor) {
	    _monitors.erase(p);
	    break;
	}
    }
    setSampledExtra();
}

void Model::addExtraNode(Node *node)
{
    if (!_is_initialized) {
	throw logic_error("Attempt to add extra node to uninitialized model");
    }
    if (node->isObserved()) {
	for (unsigned int i = 0; i < node->parents().size(); ++i) {
	    if (!node->parents()[i]->isObserved())
		throw logic_error("Cannot add observed node to initialized model");
	}
    }
    if (!node->stochasticChildren()->empty() || !node->deterministicChildren()->empty()) {
	throw logic_error("Cannot add extra node with children");
    }

    if (_extra_nodes.count(node)) {
	throw logic_error("Extra node already in model");
    }
    
    for (vector<Node const *>::const_iterator p = node->parents().begin(); 
	 p != node->parents().end(); ++p)
    {
	if (!_graph.contains(*p)) {
	throw logic_error("Extra node has parents not in model");
	}
    }

    if (!_graph.contains(node)) {
	_graph.add(node);
    }

    _extra_nodes.insert(node);
    if (_data_gen) {
	//Extra nodes are automatically sampled
	_sampled_extra.push_back(node);
    }
}


/* 
   We use construct-on-first-use for the factory lists used by model
   objects. By dynamically allocating a list, we ensure that its
   destructor is never called - the memory is simply returned to the
   OS on exit.

   This fixes a nasty exit bug.  We cannot guarantee the order that
   static destructors are called in.  Therefore, a segfault can occur
   if a module tries to remove entries from a list that has already
   been destroyed.

   See also Compiler.cc, where the same technique is used for 
   lookup tables used by the compiler.
*/

list<pair<SamplerFactory *, bool> > &Model::samplerFactories()
{
    static list<pair<SamplerFactory *, bool> > *_samplerfac =
	new list<pair<SamplerFactory *, bool> >();
    return *_samplerfac;
}

list<pair<RNGFactory *, bool> > &Model::rngFactories()
{
    static list<pair<RNGFactory *, bool> > *_rngfac =
	new list<pair<RNGFactory *, bool> >();
    return *_rngfac;
}

list<pair<MonitorFactory *, bool> > &Model::monitorFactories()
{
    static list<pair<MonitorFactory *, bool> > *_monitorfac =
	new list<pair<MonitorFactory*,bool> >();
    return *_monitorfac;
}

unsigned int Model::nchain() const
{
  return _nchain;
}

RNG *Model::rng(unsigned int chain) const
{
  return _rng[chain];
}

bool Model::setRNG(string const &name, unsigned int chain)
{
  if (chain >= _nchain)
     throw logic_error("Invalid chain number in Model::setRNG");

  list<pair<RNGFactory*, bool> >::const_iterator p;
  for (p = rngFactories().begin(); p != rngFactories().end(); ++p) {
      if (p->second) {
	  RNG *rng = p->first->makeRNG(name);
	  if (rng) {
	      _rng[chain] = rng;
	      return true;
	  }
      }
  }

  return false;
}

bool Model::setRNG(RNG *rng, unsigned int chain)
{
  if (chain >= _nchain)
     throw logic_error("Invalid chain number in Model::setRNG");

  _rng[chain] = rng;
  return true;
}

list<MonitorControl> const &Model::monitors() const
{
  return _monitors;
}

/*
bool Model::setDefaultMonitors(string const &type, unsigned int thin)
{
    list<MonitorFactory*> const &faclist = monitorFactories();

    for(list<MonitorFactory*>::const_iterator j = faclist.begin();
	j != faclist.end(); ++j)
    {
	vector <Node const *> default_nodes = (*j)->defaultNodes(this, type);
	if (!default_nodes.empty()) {
	    unsigned int start = iteration() + 1;
	    for (unsigned int i = 0; i < default_nodes.size(); ++i) {
		Monitor *monitor = (*j)->getMonitor(default_nodes[i], this,
						    type);
		if (!monitor) {
		    throw logic_error("Invalid default monitor");
		}
		addMonitor(monitor, thin);
		// Model takes ownership of default monitors
		_default_monitors.push_back(monitor);
	    }
	    return true;
	}
    }
    return false;
}

void Model::clearDefaultMonitors(string const &type)
{
    list<Monitor*> dmonitors = _default_monitors;
    for (list<Monitor*>::const_iterator p = dmonitors.begin();
	 p != dmonitors.end(); ++p) 
    {
	Monitor *monitor = *p;
	if (monitor->type() == type) {
	    _default_monitors.remove(monitor);
	    removeMonitor(monitor);
	    delete monitor;
	}
    }
    setSampledExtra();
}
*/

void Model::addNode(StochasticNode *node)
{
    _graph.add(node);
    _stochastic_nodes.push_back(node);
}

void Model::addNode(DeterministicNode *node)
{
    _graph.add(node);
}

void Model::addNode(ConstantNode *node)
{
    _graph.add(node);
}

vector<StochasticNode*> const &Model::stochasticNodes() const
{
    return _stochastic_nodes;
}
