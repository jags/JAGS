#include <config.h>
#include <Console.h>
#include <compiler/Compiler.h>
#include <compiler/parser_extra.h>
#include <compiler/ParseTree.h>
#include <model/TraceMonitor.h>
#include <model/BUGSModel.h>
#include <graph/NodeError.h>
#include <sampler/SamplerFactory.h>
#include <graph/Node.h>
#include <sarray/Range.h>
#include <sarray/SArray.h>
#include <rng/RNG.h>

#include <map>
#include <list>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <algorithm>

using std::ostream;
using std::ofstream;
using std::endl;
using std::flush;
using std::vector;
using std::string;
using std::map;
using std::list;
using std::set;
using std::pair;

/* Helper functions for dumpState */

static bool isData(Node const *node)
{
  return node->isVariable() && node->isObserved();
}

static bool isParameter(Node const *node)
{
  return node->isVariable() && !node->isObserved();
}

static bool alwaysTrue(Node const *node)
{
  return true;
}


Console::Console(ostream &out, ostream &err)
  : _out(out), _err(err), _model(0), _pdata(0), _prelations(0), 
    _pvariables(0)
{
}

Console::~Console()
{
  delete _model;
  delete _pdata;
  delete _prelations;
  if (_pvariables) {
      for (unsigned int i = 0; i < _pvariables->size(); ++i) {
	  delete (*_pvariables)[i];
      }
      delete _pvariables;
  }
}

bool Console::checkModel(FILE *file)
{
    if (_model) {
	_out << "Replacing existing model" << endl;
	clearModel();
    }
    _model = 0;
    int lineno =  parse_bugs(file, &_pvariables, &_pdata, &_prelations);
    if (lineno != 0) {
	_err << endl << "Parse error on line " << lineno << endl;
	clearModel();
	return false;
    }
    return true;
}


bool Console::compile(map<string, SArray> &data_table, unsigned int nchain,
                      bool gendata)
{
  if (nchain == 0) {
    _err << "You must have at least one chain" << endl;
    return false;
  }
  
  RNG *datagen_rng = 0;
  if (_pdata && gendata) {
    BUGSModel datagen_model(1);
    Compiler compiler(datagen_model.graph(), datagen_model.symtab(), 
		      data_table);
    _out << "Compiling data graph" << endl;
    try {
      if (_pvariables) {
	_out << "   Declaring variables" << endl;
	compiler.declareVariables(*_pvariables);
      }
      _out << "   Resolving undeclared variables" << endl;
      compiler.undeclaredVariables(_pdata);
      _out << "   Allocating nodes" << endl;
      compiler.writeRelations(_pdata);
      _out << "   Checking graph" << endl;
      datagen_model.checkGraph();
      
      vector<Node*> nodes;
      datagen_model.graph().getNodes(nodes);
      for (unsigned int i = 0; i < nodes.size(); ++i) {
	if (nodes[i]->isObserved()) {
	  vector<Node*> const &parents = nodes[i]->parents();
	  for (vector<Node*>::const_iterator p = parents.begin();
	       p != parents.end(); ++p)
	    {
	      if (!((*p)->isObserved())) {
		_err << "Invalid data graph: observed node " 
		     << nodes[i]->name(datagen_model.symtab()) 
		     << " has unobserved parent " 
		     << (*p)->name(datagen_model.symtab()) 
		     << "\n";
		return false;
	      }
	    }
	}
      }
      _out << "   Initializing" << endl;
      datagen_model.initialize();
      datagen_rng = datagen_model.rng(0);
      vector <Node*> data_nodes;
      datagen_model.graph().getSortedNodes(data_nodes);
      for (unsigned long i = 0; i < data_nodes.size(); ++i) {
	if (!data_nodes[i]->isObserved()) {
	  data_nodes[i]->randomSample(datagen_rng, 0);
	}
      }
      _out << "   Reading data back into data table" << endl;
      datagen_model.symtab().readValues(data_table, 0, alwaysTrue);
    }
    catch (NodeError except) {
      _err << "Error in node " << 
	datagen_model.symtab().getName(except.node) << "\n";
      _err << except.what() << endl;
      return false;
    }
    catch (std::runtime_error except) {
      _err << "RUNTIME ERROR:\n";
      _err << except.what() << endl;
      return false;
    }
    catch (std::logic_error except) {
      _err << "LOGIC ERROR:\n";
      _err << except.what() << '\n';
      _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
      return false;
    }
  }

  _model = new BUGSModel(nchain);
  Compiler compiler(_model->graph(), _model->symtab(), data_table);

  _out << "Compiling model graph" << endl;
  try {
    if (_pvariables) {
      _out << "   Declaring variables" << endl;
      compiler.declareVariables(*_pvariables);
    }
    if (_prelations) {
       _out << "   Resolving undeclared variables" << endl;
       compiler.undeclaredVariables(_prelations);
       _out << "   Allocating nodes" << endl;
       compiler.writeRelations(_prelations);
     }
     else {
       _err << "Nothing to compile" << endl;
       return false;
     }
     if (_model) {
       _out << "   Checking graph" << endl;
       _model->checkGraph();
       _out << "   Graph Size: " << _model->graph().size() << endl;
       if (datagen_rng) {
          // Reuse the data-generation RNG, if there is one, for chain 0 
          _model->setRNG(datagen_rng, 0);
       }
     }
     else {
       _err << "No model" << endl;
       return false;
     }
  }
  catch (NodeError except) {
    _err << "Error in node " << _model->symtab().getName(except.node) << '\n';
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    clearModel();
    return false;
  }
  return true;
}

bool Console::initialize()
{
    if (_model == 0) {
	_err << "Can't initialize. No model!" << endl;
	return false;
    }
    if (_model->graph().size() == 0) {
	_err << "Can't initialize. No nodes in graph (Have you compiled the model?)" 	 << endl;
	return false;
    }
    try {
	_model->initialize();
	_model->chooseSamplers();
	_model->addDevianceNode();
    }
    catch (NodeError except) {
	_err << "Error in node " << _model->symtab().getName(except.node) << endl;
	_err << except.what() << endl;
	return false;
    }
    catch (std::runtime_error except) {
	_err << "RUNTIME ERROR" << endl;
	_err << except.what() << endl;
	return false;
    }
    catch (std::logic_error except) {
	_err << "LOGIC ERROR" << endl;
	_err << except.what() << endl;
	_err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
	return false;
    }
    return true;
}

bool Console::setParameters(map<string, SArray> const &init_table,
			    unsigned int chain)
{
  if (_model == 0) {
    _err << "Can't set initial values. No model!" << endl;    
    return false;
  }
  if (chain == 0 || chain > nchain()) {
    _err << "Invalid chain number" << endl;
    return false;
  }

  try {
    _model->setParameters(init_table, chain - 1);
  }
  catch (NodeError except) {
    _err << "Error in node " << _model->symtab().getName(except.node) << '\n';
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    clearModel();
    return false;
  }
  return true;
}

bool Console::setRNGname(string const &name, unsigned int chain)
{
  if (_model == 0) {
    _err << "Can't set RNG name. No model!" << endl;    
    return false;
  }
  try {
    bool name_ok = _model->setRNG(name, chain - 1);
    if (!name_ok) {
       _err << "WARNING: RNG name " << name << " not found\n";
    }
  }
  catch (NodeError except) {
    _err << "Error in node " << _model->symtab().getName(except.node) << '\n';
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    clearModel();
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    clearModel();
    return false;
  }
  return true;
}

bool Console::update(unsigned long n)
{
    if (_model == 0) {
	_err << "Can't update. No model!" << endl;    
	return false;
    }
    if (!_model->isInitialized()) {
	_err << "Model not initialized" << endl;
	return false;
    }
    if (!_model->canSample()) {
	_err << "Model has no samplers" << endl;
        return false;
    }
    try {
	_model->update(n);
    }
    catch (NodeError except) {
	_err << "Error in node " << _model->symtab().getName(except.node) << '\n';
	_err << except.what() << endl;
	return false;
    }
    catch (std::runtime_error except) {
	_err << "RUNTIME ERROR:\n";
	_err << except.what() << endl;
	return false;
    }
    catch (std::logic_error except) {
	_err << "LOGIC ERROR:\n";
	_err << except.what() << '\n';
	_err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
	return false;
    }
    return true;
}

unsigned long Console::iter(unsigned int chain) const
{
  if (!_model) {
    return 0;
  }
  else {
    return _model->iteration(chain - 1);
  }
}

bool Console::setMonitor(string const &name, Range const &range,
			 unsigned int thin)
{
  if (!_model) {
    _out << "Can't set monitor. No model!" << endl;    
    return false;
  }
  
  try {
    if (isNULL(range)) {
      NodeArray *array = _model->symtab().getVariable(name);
      if (array) {
	return setMonitor(name, array->range(), thin);
      }
      else {
	_err << "Variable " << name << " not found" << endl;
	return false;
      }
    }
    else {
      string message;
      bool status = _model->setMonitor(name, range, thin, message);
      if (!status) {
         _err << message << endl;
         return false;
      }
    }
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    return false;
  }
  return true;
}

bool Console::clearMonitor(string const &name, Range const &range)
{
  if (!_model) {
    _out << "Can't clear monitor. No model!" << endl;    
    return false;
  }

  try {
    Node *node = 0;
    string message;
    if (isNULL(range)) {
      NodeArray *array = _model->symtab().getVariable(name);
      if (array) {
        node = _model->getNode(name, array->range(), message);
      }
      else {
	_err << "Variable " << name << " not found" << endl;
	return false;
      }
    }
    else {
      Node *node = _model->getNode(name, range, message);
    }

    if (node) {
      _model->clearMonitor(node);	  
    }
    else {
      _err << message << endl;
      return false;
    }
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:" << endl;
    _err << except.what() << endl;
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:" << endl;
    _err << except.what() << endl;
    return false;
  }
  return true;
}

void Console::clearModel()
{
    /*
    if (_pvariables) {
	for (unsigned int i = 0; i < _pvariables->size(); ++i) {
	    delete (*_pvariables)[i];
	}
	delete _pvariables;
	_pvariables = 0;
    }
    delete _pdata; _pdata = 0;
    delete _prelations; _prelations = 0;
    */
    _out << "Deleting model" << endl;
    delete _model; 
    _model = 0;
}

bool Console::dumpState(map<string,SArray> &data_table,
			string &rng_name,
			DumpType type, unsigned int chain)
{
  if (!_model) {
    _err << "No model" << endl;
    return false;
  }
  if (_model->symtab().size() == 0) {
    _err << "Symbol table is empty" << endl;
    return false;
  }
  if (chain == 0 || chain > nchain()) {
    _err << "Invalid chain number" << endl;
    return false;
  }
  bool (*selection)(Node const *);
  switch (type) {
  case DUMP_PARAMETERS:
    selection = isParameter;
    break;
  case DUMP_DATA:
    selection = isData;
    break;
  case DUMP_ALL:
    selection = alwaysTrue;
    break;
  }

  try {
    _model->symtab().readValues(data_table, chain - 1, selection);
    
    if (type == DUMP_PARAMETERS || type == DUMP_ALL) {
      
      vector<int> rngstate;
      if (_model->rng(chain - 1)) {
	_model->rng(chain - 1)->getState(rngstate);
	
	Index dimrng(1);
	dimrng[0] = rngstate.size();
	SArray rngsarray(dimrng);
	rngsarray.setDiscreteValued(true);
	for (unsigned i = 0; i < rngstate.size(); ++i) {
	  rngsarray.setValue(rngstate[i], i);
	}
	
	data_table.insert(pair<string, SArray>(".RNG.state",rngsarray));
	rng_name = _model->rng(chain - 1)->name();
      }
    }
  }
  catch (NodeError except) {
    _err << "Error in node " <<
      _model->symtab().getName(except.node) << "\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    return false;
  }
  
  return true;
}

bool Console::coda(ofstream &index, vector<ofstream*> &output)
{
  if (!_model) {
     _err << "Can't dump CODA output. No model!" << endl;
     return false;
  }

  vector<Node const *> nodes;
  list<TraceMonitor*> const &monitors = _model->monitors(0);
  for (list<TraceMonitor*>::const_iterator i = monitors.begin();
       i != monitors.end(); i++) 
    {
      nodes.push_back((*i)->node());
    }

  try {
    _model->coda(nodes, index, output);
  }
  catch (NodeError except) {
    _err << "Error in node " <<
      _model->symtab().getName(except.node) << "\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:\n";
      _err << except.what() << '\n';
      _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;            return false;
  }
  
  return true;
}

bool Console::getMonitoredValues(map<string,SArray> &data_table,
				 unsigned int chain)
{
  if(chain == 0 || chain > nchain()) {
    _err << "Invalid chain number" << endl;
    return false;
  }

  try {
    list<TraceMonitor*> const &monitors = _model->monitors(chain - 1);
    list<TraceMonitor*>::const_iterator p;
    for (p = monitors.begin(); p != monitors.end(); ++p) {
      TraceMonitor const *monitor = *p;
      if (monitor->size() > 0) {
	Node const *node = monitor->node();
	string name = _model->symtab().getName(node);

	//The new SArray has the same dimensions as the
	//monitored node, plus an extra one for the 
	//iterations. We put the extra dimension first.
	unsigned int niter = monitor->size();
	unsigned int ndim = node->dim(false).size();
	Index dim(ndim + 1);
	dim[0] = niter;
	unsigned long length = 1;
	for (unsigned int i = 1; i <= ndim; ++i) {
	    dim[i] = node->dim(false)[i-1];
	    length *= dim[i];
	}
	//Create a new SArray and insert it into the table
	SArray ans(dim);
	double *values = new double[length * niter];
	double const *monitor_values = monitor->values();
	for (unsigned int i = 0; i < length; ++i) {
	    for (unsigned int j = 0; j < niter; ++j) {
		values[niter * i + j] = monitor_values[length * j + i];
	    }
	}
	ans.setValue(values, length*niter);
	delete [] values;
	data_table.insert(pair<string,SArray>(name, ans));
      }
    }
  }
  catch (NodeError except) {
    _err << "Error in node " <<
      _model->symtab().getName(except.node) << "\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::runtime_error except) {
    _err << "RUNTIME ERROR:\n";
    _err << except.what() << endl;
    return false;
  }
  catch (std::logic_error except) {
    _err << "LOGIC ERROR:\n";
    _err << except.what() << '\n';
    _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;
    return false;
  }
}

bool Console::coda(string const &name, Range const &range,
		   ofstream &index, vector<ofstream*> &output)
{
  if (!_model) {
     _err << "Can't dump CODA output. No model!" << endl;
     return false;
  }

  Range frange = range;
  NodeArray *array = _model->symtab().getVariable(name);
  if (array) {
    if (isNULL(frange)) {
      frange = array->range();
    }
    if (array->range().contains(frange)) {
      Node *node = array->getSubset(frange);
      if (node) {
	try {
	  _model->coda(vector<Node const*>(1, node), index, output);
	}
	catch (NodeError except) {
	  _err << "Error in node " <<
	    _model->symtab().getName(except.node) << "\n";
	  _err << except.what() << endl;
	  return false;
	}
	catch (std::runtime_error except) {
	  _err << "RUNTIME ERROR:\n";
	  _err << except.what() << endl;
	  return false;
	}
	catch (std::logic_error except) {
	  _err << "LOGIC ERROR:\n";
	  _err << except.what() << '\n';
	  _err << "Please send a bug report to " << PACKAGE_BUGREPORT << endl;            return false;
	}

	return true;
      }
      else {
	_err << "Node is not being monitored" << endl;
	return false;
      }
    }
    else {
      _err << "Requested invalid subset of node " << name << endl;
      return false;
    }
  }
  else {
    _err << name << " not found" << endl;
    return false;
  }
}

/*
bool Console::coda(string const &name, ofstream &index, ofstream &output)
{
  if (!_model) {
    _err << "Can't dump CODA output. No model!" << endl;
    return false;
  }
  
  NodeArray *array = _model->symtab().getVariable(name);
  if (array) {
    vector<Node const *> nodes(1, array->getSubset(array->range()));
    _model->coda(nodes, index, output);
    return true;
  }
  else {
    _err << name << " not found" << endl;
    return false;
  }
}
*/
  
BUGSModel const *Console::model()
{
  return _model;
}

/*
void jags_set_seed(unsigned int seed)
{
  // Given a single seed, calculates seeds for the Marsaglia
  // Multicarry random number generator
  unsigned int mmseeds[2];
  for(int j = 0; j < 2; j++) {
    seed = (69069 * seed + 1);
    mmseeds[j] = seed;
  }
  // Fix up zero seeds 
  if (mmseeds[0] == 0)
    mmseeds[0] = 1;
  if (mmseeds[1] == 0)
    mmseeds[1] = 1;

  set_seed(mmseeds[0], mmseeds[1]);
}
*/

unsigned int Console::nchain() const
{
  if (_model == 0) {
    return 0;
  }
  else {
    return _model->nchain();
  }
}

/*
bool Console::setRNGseed(SArray const &seed, unsigned int chain)
{
  if (_model == 0) {
    _err << "Can't set RNG seed. No model!" << endl;    
    return false;
  }
  if (seed.length() != 1) {
    _err << "Seed must be a single integer\n";
    return false;
  }
  if (chain == 0 || chain > _model->nchain()) {
    _err << "Invalid chain number\n";
    return false;
  }

  int iseed = static_cast<int>(*seed.value());
  _model->rng(chain-1)->init(iseed);
  return true;
}

bool Console::setRNGstate(SArray const &state, unsigned int chain)
{
  if (_model == 0) {
    _err << "Can't set RNG state. No model!" << endl;    
    return false;
  }
  vector<int>(istate);
  double const *value = state.value();
  for (unsigned int i = 0; i < state.length(); ++i) {
    istate.push_back(static_cast<int>(value[i]));
  }

  if (_model->rng(chain-1)->setState(istate) == false) {
    _err << "Invalid RNG state\n";
    return false;
  }
  return true;
}
*/
