

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
    
    Graph full_graph;
    for (unsigned int i = 0; i < _nodes.size(); ++i) {
	full_graph.insert(_nodes[i]);
    }

    GraphMarks marks(full_graph);
    Graph sample_graph;

    // Add observed stochastic nodes to the sample graph and mark
    // the informative nodes

    vector<Node const*> informative;
    vector<StochasticNode*>::const_iterator p;
    for (p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	if (isObserved(*p)) {
	    sample_graph.insert(*p);
	    informative.push_back(*p);
	}
    }
    marks.markAncestors(informative, 1);

    for (p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	if (isObserved(*p)) {
	    marks.mark(*p, 2);
	}
    }

    //Triage on marked nodes. We do this twice: once for stochastic
    //nodes and once for all nodes.

    list<StochasticNode*> slist; //List of nodes to be sampled
    for(p = _stochastic_nodes.begin(); p != _stochastic_nodes.end(); ++p) {
	if (marks.mark(*p) == 1) {
	    //Unobserved stochastic nodes: to be sampled
	    slist.push_back(*p); 
	}
    }

    for (vector<Node*>::const_iterator j = _nodes.begin();
	 j != _nodes.end(); ++j) 
    {
	switch(marks.mark(*j)) {
	case 0:
	    _extra_nodes.push_back(*j);
	    break;
	case 1: case 2:
	    sample_graph.insert(*j);
	    break;
	}
    }

    // Traverse the list of samplers, selecting nodes that can be sampled
    list<pair<SamplerFactory *, bool> > const &sf = samplerFactories();
    for(list<pair<SamplerFactory *, bool> >::const_iterator q = sf.begin();
	q != sf.end(); ++q) 
    {
	if (!q->second) continue;

	vector<Sampler*> svec = q->first->makeSamplers(slist, sample_graph);
	while (!svec.empty()) {
	    for (unsigned int i = 0; i < svec.size(); ++i) {

		vector<StochasticNode*> const &nodes = svec[i]->nodes();
		for (unsigned int j = 0; j < nodes.size(); ++j) {
		    /* FIXME: This is a potential bottleneck if slist
		       is large */
		    list<StochasticNode*>::iterator p = 
			find(slist.begin(), slist.end(), nodes[j]);
		    if (p == slist.end()) {
			throw logic_error("Unable to find sampled node");
		    }
		    slist.erase(p);
		}
		_samplers.push_back(svec[i]);
	    }
	    svec = q->first->makeSamplers(slist, sample_graph);
	}
    }
  
    // Make sure we found a sampler for all the nodes
    if (!slist.empty()) {
	throw NodeError(*slist.begin(),
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
    for (vector<Node *>::const_iterator p = _extra_nodes.begin(); 
	 p != _extra_nodes.end(); ++p)
    {
	egraph.insert(*p);
    }
    //Mark the ancestors of all monitored nodes in this graph
    GraphMarks emarks(egraph);
    vector<Node const*> monitored_nodes;
    for (list<MonitorControl>::const_iterator p = _monitors.begin();
	 p != _monitors.end(); ++p)
    {
	vector<Node const*> const &pnodes = p->monitor()->nodes();
	for (vector<Node const*>::const_iterator i = pnodes.begin();
	     i != pnodes.end(); ++i)
	{
	    if (egraph.contains(*i)) {
		emarks.mark(*i, 1);
		monitored_nodes.push_back(*i);
	    }
	}
    }
    emarks.markAncestors(monitored_nodes, 1);

    //Add marked nodes to the vector of sampled extra nodes
    _sampled_extra.clear();
    for(vector<Node *>::const_iterator p = _extra_nodes.begin();
	p != _extra_nodes.end(); ++p)
    {
	if (emarks.mark(*p)) {
	    _sampled_extra.push_back(*p);
	}
    }
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
	      /* NO! RNGs are owned by the factory, not the model
	      if (_rng[chain])
		  delete _rng[chain];
	      */
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

void Model::addNode(StochasticNode *node)
{
    _nodes.push_back(node);
    _stochastic_nodes.push_back(node);
}

void Model::addNode(DeterministicNode *node)
{
    _nodes.push_back(node);
}

void Model::addNode(ConstantNode *node)
{
    _nodes.push_back(node);
}

vector<StochasticNode*> const &Model::stochasticNodes() const
{
    return _stochastic_nodes;
}

    vector<Node*> const &Model::nodes() const
    {
	return _nodes;
    }

} //namespace jags
