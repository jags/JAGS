#include <config.h>
#include <model/BUGSModel.h>
#include <model/Monitor.h>
#include <model/NodeArray.h>
#include <model/MonitorFactory.h>
#include <graph/DevianceNode.h>
#include <graph/StochasticNode.h>
#include <graph/GraphMarks.h>
#include <graph/Node.h>
#include <rng/RNG.h>
#include <sampler/Sampler.h>

#include <list>
#include <utility>
#include <stdexcept>
#include <fstream>
#include <cmath>

using std::vector;
using std::ofstream;
using std::list;
using std::pair;
using std::string;
using std::logic_error;
using std::runtime_error;

/* 
   Nodes accessible to the user in a BUGSModel are identified
   by a variable name and range of indices 
*/
typedef pair<string, Range> NodeId;

BUGSModel::BUGSModel(unsigned int nchain)
  : Model(nchain), _symtab(graph(), nchain)
{
}

BUGSModel::~BUGSModel()
{
    for (list<Monitor*>::iterator i = _bugs_monitors.begin();
	 i != _bugs_monitors.end(); ++i)
    {
	Monitor *monitor = *i;
	/*
	removeMonitor(monitor);
	_bugs_monitors.erase(i);
	*/
	delete monitor;
    }

}

SymTab &BUGSModel::symtab()
{
    return _symtab;
}

Node *BUGSModel::getNode(string const &name, Range const &target_range,
			 string &message)
{
    NodeArray *array = _symtab.getVariable(name);
    if (!array) {
	message = string("Unknown variable ") + name;
	return 0;
    }

    Range range = target_range;
    if (isNULL(range)) {
	range = array->range();
    }
    else if (!array->range().contains(target_range)) {
	//FIXME. What if target_range has wrong dimension?
	message = string("Invalid range ") + name + print(target_range);
	return 0;
    }
    message.clear();
    unsigned int NNode = graph().size();
    Node *node = array->getSubset(range, graph());
    if (graph().size() != NNode) {
	addExtraNode(node); // Node was newly allocated
    }
    return node;
}

//FIXME Should be a util
#include <util/nainf.h>
static void writeDouble(double x, std::ostream &out)
{
    if (x == JAGS_NA) {
	out << "NA";
    }
    else if (jags_isnan(x)) {
	out << "NaN";
    }
    else if (!jags_finite(x)) {
	if (x > 0) 
	    out << "Inf";
	else 
	    out << "-Inf";
    }
    else {
	out << x;
    }
}

static void CODA(vector<pair<NodeId, Monitor const*> > const &monitors,
		 ofstream &index,  vector<ofstream*> const &output)
{
    //FIXME: Check streams are open and set pointer to beginning.

    unsigned int lineno = 0;
    vector<pair<NodeId, Monitor const*> >::const_iterator p;
    for (p = monitors.begin(); p != monitors.end(); ++p) {
	
	Monitor const *monitor = p->second;
	if (monitor->type() == "trace") {
	    
	    string const &name = p->first.first;
	    Range const &uni_range = p->first.second;
	    Range multi_range = uni_range;
	    if (isNULL(multi_range)) {
		multi_range = Range(monitor->node()->dim());
	    }
	    else if (multi_range.dim(true) != monitor->node()->dim()) {
		throw logic_error("Range does not match Node dimension in CODA");
	    }

	    unsigned int nvar = monitor->node()->length();
	    if (nvar != 1) {
		// Multivariate node
		for (unsigned int offset = 0; offset < nvar; ++offset) 
		{
		    index << name << print(multi_range.leftIndex(offset))
			  << " "  << lineno + 1 << "  "  
			  << lineno + monitor->niter() << '\n';
		    lineno += monitor->niter();
		}
	    }
	    else {
		// Univariate node 
		index << name << print(uni_range) << "  " 
		      << lineno + 1 << "  " 
		      << lineno + monitor->niter() << '\n';
		lineno += monitor->niter();
	    }
	    //Write output files
	    for (unsigned int ch = 0; ch < output.size(); ++ch) {
		ofstream &out = *output[ch];
		vector<double> const &y = monitor->value(ch);
		for (unsigned int offset = 0; offset < nvar; ++offset) {
		    unsigned int iter = monitor->start();
		    for (unsigned int k = 0; k < monitor->niter(); k++) {
			out << iter << "  ";
			writeDouble(y[k * nvar + offset], out);
			out << '\n';
			iter += monitor->thin();
		    }
		}
	    }
	}
    }
}

void BUGSModel::coda(vector<NodeId> const &nodes, ofstream &index, 
		     vector<ofstream*> const &output, string &warn)
{
    if (output.size() != nchain()) {
        throw logic_error("Wrong number of output streams in BUGSModel::coda");
    }

    warn.clear();

    // Create a new vector with only the selected nodes in it 
    vector<pair<NodeId, Monitor const*> > dump_nodes;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	string msg;
	string const &name = nodes[i].first;
	Range const &range = nodes[i].second;
	Node *node = getNode(name, range, msg);
	if (!node) {
	    warn.append(msg);
	}
	else {
	    list<Monitor*>::const_iterator p = _bugs_monitors.begin();
	    for (; p != _bugs_monitors.end(); ++p) {
		if ((*p)->node() == node && (*p)->type() == "trace") {
		    pair<NodeId, Monitor const*> newpair(nodes[i], *p);
		    dump_nodes.push_back(newpair);		    
		    break;
		}
	    }
	    if (p == _bugs_monitors.end()) {
		string msg  = string("No trace monitor for node ") + name + 
		    print(range);
		warn.append(msg);
	    }
	}
    }

    if (dump_nodes.empty()) {
	warn.append("There are no monitored nodes\n");
    }
    
    CODA(dump_nodes, index, output);
}

void BUGSModel::coda(ofstream &index, vector<ofstream*> const &output,
                     string &warn)
{
    if (output.size() != nchain()) {
        throw logic_error("Wrong number of output streams in BUGSModel::coda");
    }
    
    warn.clear();
    
    vector<pair<NodeId,Monitor const*> > dump_nodes;
    for(list<Monitor*>::const_iterator p = _bugs_monitors.begin();
	p != _bugs_monitors.end(); ++p)
    {
	Monitor const *monitor = *p;
	if (monitor->type() == "trace") {
	    NodeId id = _node_map.find(monitor->node())->second;
	    dump_nodes.push_back(pair<NodeId,Monitor const*>(id,monitor));
	}
    }
    
    if (dump_nodes.empty()) {
	warn.append("There are no monitored nodes\n");
    }
 
    CODA(dump_nodes, index, output);
}

void BUGSModel::addDevianceNode()
{
    NodeArray const *array = _symtab.getVariable("deviance");
    if (array)
	return; //Deviance already defined by user

    _symtab.addVariable("deviance", vector<unsigned int>(1,1));
    NodeArray *deviance = _symtab.getVariable("deviance");
    vector<Node*> nodes;
    graph().getNodes(nodes);
    std::set<StochasticNode const *> parameters;
    for (vector<Node*>::const_iterator p = nodes.begin(); p != nodes.end(); ++p)
    {
	if ((*p)->isObserved()) {
	    StochasticNode *snode = dynamic_cast<StochasticNode*>(*p);
	    if (snode)
		parameters.insert(snode);
	}
    }
    if (!parameters.empty()) {
	//Can't construct a deviance node with no parameters
	DevianceNode *dnode = new DevianceNode(parameters);
	addExtraNode(dnode);
	deviance->insert(dnode, vector<unsigned int>(1,1));
    }
}

void BUGSModel::setParameters(std::map<std::string, SArray> const &param_table,
			      unsigned int chain)
{
    _symtab.writeValues(param_table, chain);


    //Strip off .RNG.seed (user-supplied random seed)
    if (param_table.find(".RNG.seed") != param_table.end()) {
	if (rng(chain) == 0) {
	    throw runtime_error(".RNG.seed supplied but RNG type not set");
	}
	SArray const &seed = param_table.find(".RNG.seed")->second;
	if (seed.length() != 1) {
	    throw runtime_error(".RNG.seed must be a single integer");
	}
	if (seed.value()[0] < 0) {
	    throw runtime_error(".RNG.seed must be non-negative");
	}
	int iseed = static_cast<int>(seed.value()[0]);
	rng(chain)->init(iseed);
    }

    //Strip off .RNG.state (saved state from previous run)
    if (param_table.find(".RNG.state") != param_table.end()) {
	if (rng(chain) == 0) {
	    throw runtime_error(".RNG.state supplied, but RNG type not set");
	}
	SArray const &state = param_table.find(".RNG.state")->second;
	vector<int>(istate);
	//double const *value = state.value();
	vector<double> const &value = state.value();
	for (unsigned int i = 0; i < state.length(); ++i) {
	    istate.push_back(static_cast<int>(value[i]));
	}
	if (rng(chain)->setState(istate) == false) {
	    throw runtime_error("Invalid .RNG.state");
	}
    }
}


bool BUGSModel::setMonitor(string const &name, Range const &range,
			   unsigned int thin, string const &type)
{
    string msg;
    Node *node = getNode(name, range, msg);
    if (!node) {
	return false;
    }

    for (list<Monitor*>::const_iterator i = _bugs_monitors.begin();
	 i != _bugs_monitors.end(); ++i)
    {
	if ((*i)->node() == node && (*i)->type() == type)
	    return false; //Node is already being monitored.
    }

    Monitor *monitor = 0;

    list<MonitorFactory*> const &faclist = monitorFactories();
    for(list<MonitorFactory*>::const_iterator j = faclist.begin();
	j != faclist.end(); ++j)
    {
	monitor = (*j)->getMonitor(node, this, iteration() + 1, thin, type);
	if (monitor)
	    break;
    }

    if (monitor) {
	addMonitor(monitor);
	_bugs_monitors.push_back(monitor);
	_node_map.insert(pair<Node const*,NodeId>(node, NodeId(name,range)));
	return true;
    }
    else {
	return false;
    }
}

bool BUGSModel::deleteMonitor(string const &name, Range const &range,
			      string const &type)
{
    string msg;
    Node *node = getNode(name, range, msg);
    if (!node) {
	return false;
    }

    for (list<Monitor*>::iterator i = _bugs_monitors.begin();
	 i != _bugs_monitors.end(); ++i)
    {
	if ((*i)->node() == node && (*i)->type() == type) {
	    Monitor *monitor = *i;
	    removeMonitor(monitor);
	    _bugs_monitors.erase(i);
	    delete monitor;
	    return true;
	}
    }
    return false;
}

void BUGSModel::samplerNames(vector<vector<string> > &sampler_names) const
{
    sampler_names.clear();
    sampler_names.reserve(_samplers.size());

    for (unsigned int i = 0; i < _samplers.size(); ++i) {

	vector<string> names;	
	vector<StochasticNode *> const &nodes = _samplers[i]->nodes();
	names.reserve(nodes.size()+1);
	
	names.push_back(_samplers[i]->name());
	for (unsigned int j = 0; j < nodes.size(); ++j) {
	    names.push_back(_symtab.getName(nodes[j]));
	}
	sampler_names.push_back(names);
    }    
}
