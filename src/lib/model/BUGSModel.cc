#include <config.h>
#include <model/BUGSModel.h>
#include <model/TraceMonitor.h>
#include <model/NodeArray.h>
#include <graph/DevianceNode.h>
#include <graph/StochasticNode.h>
#include <graph/GraphMarks.h>
#include <graph/Node.h>
#include <rng/RNG.h>

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
using std::set;

BUGSModel::BUGSModel(unsigned int nchain)
  : Model(nchain), _symtab(nchain)
{
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
  else if (!array->range().contains(target_range)) {
    message = string("Invalid range ") + name + print(target_range);
    return 0;
  }
  else {
    message.clear();
    return array->getSubset(target_range);
  }
}

bool BUGSModel::setMonitor(string const &name, Range const &range, 
			   unsigned int thin, string &message)
{
  Node *node = getNode(name, range, message);
  if (node) {
    _node_map.insert(pair<Node const*, pair<string,Range> >(node, pair<string,Range>(name,range)));
    Model::setMonitor(node, thin);
    return true;
  }
  else {
    return false;
  }

}

/*
TraceMonitor const *BUGSModel::getMonitor(string const &name, Range const &range, unsigned int chain)
{
  Node *get
  NodeArray *array = _symtab.getVariable(name);
  if (!array) {
    throw runtime_error(string("Unknown variable ") + print(range));
  }
  if (array->range().contains(range)) {
    Node const *node = array->getSubset(range);
    if(node)
      Model::clearMonitor(node);
  }
  else {
    throw runtime_error(string("Requested invalid subset ") + name
			+ print(range));
  }
}
*/

static void writeDouble(double x, std::ostream &out)
{
  using namespace std;

  if (x == JAGS_NA) {
    out << "NA";
  }
  else if (isnan(x)) {
    out << "NaN";
  }
  else if (isinf(x)) {
    if (x > 0) 
      out << "Inf";
    else 
      out << "-Inf";
  }
  else {
    out << x;
  }
}

void BUGSModel::coda(vector<Node const*> const &nodes, ofstream &index, 
		     vector<ofstream*> &output)
{
    /* 
       Dump monitor for given node in CODA format.
    */

    if (output.size() != nchain()) {
        throw logic_error("Wrong number of output streams in BUGSModel::coda");
    }

    long lineno = 0;
    for (unsigned int i = 0; i < nodes.size(); i++) {

	Node const *node = nodes[i];
	string const &name = _node_map[node].first;
	Range const &range = _node_map[node].second;
	int nvar = node->length();

	list<TraceMonitor*>::const_iterator j;
	//Write index file
	for (j = monitors(0).begin(); j != monitors(0).end(); ++j) {
	    TraceMonitor const *monitor = *j;
	    if (monitor->node() == node) {
		if (nvar != 1) {
		    /* Multivariate node */
		    for (long offset = 0; offset < nvar; ++offset) {
			index << name << print(range.leftIndex(offset))
			      << " "  
			      << lineno + 1 << "  "  
			      << lineno + monitor->size() << '\n';
		        lineno += monitor->size();
		    }
		}
		else {
		    /* Univariate node */
		    index << _symtab.getName(node) << "  " 
			  << lineno + 1 << "  " 
			  << lineno + monitor->size() << '\n';
		    lineno += monitor->size();
		}
	    }
	}
	//Write output files
	for (unsigned int ch = 0; ch < nchain(); ++ch) {
	    for (j = monitors(ch).begin(); j != monitors(ch).end(); ++j) {
		TraceMonitor const *monitor = *j;
		if (monitor->node() == node) {
		    double const *y = monitor->values();
		    for (long offset = 0; offset < nvar; ++offset) {
			long iter = monitor->start();
			for (int k = 0; k < monitor->size(); k++) {
			    *(output[ch]) << iter << "  ";
			    writeDouble(y[k * nvar + offset], *output[ch]);
			    *(output[ch]) << '\n';
			    iter += monitor->thin();
			}
		    }
		}
	    }
	}
    }
}

void BUGSModel::addDevianceNode()
{
  NodeArray const *array = _symtab.getVariable("deviance");
  if (array)
    return; //Deviance already defined by user

  _symtab.addVariable("deviance", Index(1));
  NodeArray *deviance = _symtab.getVariable("deviance");
  vector<Node*> nodes;
  graph().getNodes(nodes);
  std::set<StochasticNode*> parameters;
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
     deviance->insert(dnode, Index(1));
     addExtraNode(dnode);
  }
}

static void addDescendants(Node *node, Graph &descendants)
{
   /* Adds node and its descendants to the given graph */

   if (descendants.contains(node)) {
      return;
   }

   descendants.add(node);

   set<Node*> const &children = node->children();
   set<Node*>::const_iterator p;
   for (p = children.begin(); p != children.end(); ++p) {
       addDescendants(*p, descendants);
   }
}

void BUGSModel::setParameters(std::map<std::string, SArray> const &param_table,
			      unsigned int chain)
{
    set<Node*> setnodes;
    _symtab.writeValues(param_table, chain, setnodes);


    //Strip off .RNG.seed (user-supplied random seed)
    if (param_table.find(".RNG.seed") != param_table.end()) {
      if (rng(chain) == 0) {
	throw runtime_error(".RNG.seed supplied RNG type not set");
      }
      SArray const &seed = param_table.find(".RNG.seed")->second;
      if (seed.length() != 1) {
	throw runtime_error(".RNG.seed must be a single integer");
      }
      if (*seed.value() < 0) {
	throw runtime_error(".RNG.seed must be non-negative");
      }
      int iseed = static_cast<int>(*seed.value());
      /*
      if (rng(chain) == 0) {
          getNextRNG();
      }
      */
      rng(chain)->init(iseed);
    }

    //Strip off .RNG.state (saved state from previous run)
    if (param_table.find(".RNG.state") != param_table.end()) {
      if (rng(chain) == 0) {
         throw runtime_error(".RNG.state supplied, but RNG type not set");
      }
      SArray const &state = param_table.find(".RNG.state")->second;
      vector<int>(istate);
      double const *value = state.value();
      for (unsigned int i = 0; i < state.length(); ++i) {
	istate.push_back(static_cast<int>(value[i]));
      }
      if (rng(chain)->setState(istate) == false) {
	throw runtime_error("Invalid .RNG.state");
      }
    }

    //set<Node*>::const_iterator p;

    /*
    std::cout << "Setting values\n";
    for (p = setnodes.begin(); p != setnodes.end(); p++) {
        std::cout << _symtab.getName(*p) << " = ";
        for (unsigned int i = 0; i < (*p)->length(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << (*p)->value(chain)[i];
        }
        std::cout << "\n";
    }
    */

    /*
    GraphMarks marks(graph());
    for (p = setnodes.begin(); p != setnodes.end(); p++) {
	marks.markDescendants(*p,1);
    }
    for (p = setnodes.begin(); p != setnodes.end(); p++) {
	marks.mark(*p, 0);
    }
    */

    /*
    FIXME: This is all very well in principle, but setParameters
    may be called before the model is initialized, in which
    case it may be illegal to update the children. 

    Also, we need to work out how far to take the updating.
    Is it just the immediate deterministic descendants? Or more?
    vector<Node*> nodes;
    graph().getSortedNodes(nodes);
    for (unsigned long i = 0; i < nodes.size(); ++i) {
	if (!nodes[i]->isObserved() && marks.mark(nodes[i]) == 1) {
            std::cout << "Updating\n";
            std::cout << nodes[i]->name(_symtab) << " = ";
	    nodes[i]->deterministicSample(chain);
            for (unsigned int j = 0; j < nodes[i]->length(); ++j) {
                if (j > 0) std::cout << ",";
                std::cout << nodes[i]->value(chain)[j];
            }
            std::cout << "\n";
        }
    }
    */
}

