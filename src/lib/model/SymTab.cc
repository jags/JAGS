#include <config.h>
#include <model/SymTab.h>
#include <graph/MixtureNode.h>

#include <cfloat>
#include <string>
#include <stdexcept>
#include <utility>
#include <set>

using std::vector;
using std::map;
using std::pair;
using std::string;
using std::runtime_error;
using std::logic_error;
using std::set;

SymTab::SymTab(unsigned int nchain)
  : _nchain(nchain)
{
}

void SymTab::addVariable(string const &name, Index const &dim)
{
  if (_varTable.find(name) != _varTable.end()) {
    string msg("Name ");
    msg.append(name);
    msg.append(" already in use in symbol table");
    throw runtime_error(msg);
  }

  NodeArray *array = new NodeArray(name, dim, _nchain);
  _varTable[name] = array;
}

NodeArray* SymTab::getVariable(string const &name) const
{
  map<string, NodeArray*>::const_iterator p =  _varTable.find(name);
  
  if (p == _varTable.end()) {
    return 0;
  }
  else {
    return p->second;
  }
}

void SymTab::getNodes(std::vector<Node*> &nodes)
{
  map<string, NodeArray*>::const_iterator p(_varTable.begin());

  for (; p != _varTable.end(); ++p) {
    Graph const &graph = p->second->graph();
    graph.getNodes(nodes);
  }
}

void SymTab::writeData(std::map<std::string, SArray> const &data_table)
{
  for(map<string, SArray>::const_iterator p(data_table.begin()); 
      p != data_table.end(); ++p) {
    NodeArray *array = getVariable(p->first);
    if (array) {
      if (array->range().dim(false) != p->second.dim(false)) {
	string msg("Dimension mismatch in values supplied for ");
	msg.append(p->first);
	throw runtime_error(msg);
      }
      array->setData(p->second);
    }
  }
}


void SymTab::writeValues(std::map<std::string, SArray> const &data_table,
		         unsigned int chain, set<Node*> &setnodes)
{
    for(map<string, SArray>::const_iterator p(data_table.begin()); 
	p != data_table.end(); ++p) {
        set<Node*> psetnodes;
	NodeArray *array = getVariable(p->first);
	if (array) {
	    if (array->range().dim(false) != p->second.dim(false)) {
		string msg("Dimension mismatch in values supplied for ");
		msg.append(p->first);
		throw runtime_error(msg);
	    }
	    array->setValue(p->second, chain, psetnodes);
            setnodes.insert(psetnodes.begin(), psetnodes.end());
	}
    }
}

void SymTab::readValues(map<string, SArray> &data_table, 
		        unsigned int chain,
                        bool (*condition)(Node const *)) const
{
  if (chain > _nchain) 
    throw logic_error("Invalid chain in SymTab::readValues");

  map<string, NodeArray*>::const_iterator p;
  for (p = _varTable.begin(); p != _varTable.end(); ++p) {
    /* Create a new SArray to hold the values from the symbol table */
    SArray read_values(p->second->range().dim(false));
    p->second->getValue(read_values, chain, condition);
    /* Only write to the data table if we can find at least one
       non-missing value */
    for (unsigned long i = 0; i < read_values.length(); ++i) {
      if (read_values.value()[i] != JAGS_NA) {
	string const &name = p->first;
	if (data_table.find(name) != data_table.end()) {
	  data_table.erase(name);
	}
	data_table.insert(pair<string,SArray>(name, read_values));
      }
    }
  }
}

unsigned int SymTab::size() const
{
  return _varTable.size();
}

void SymTab::clear()
{
  _varTable.clear();
}

static vector<string> cutBUGSSubsetName(string const &name)
{
  /* 
     Takes the BUGS language name of an array subset, e.g. "foo[a,b,c]",
     and cuts it up into substrings ("foo" "[" "a" "," "b" "," "c" "]")
  */
  vector<string> cut;
  char const *c = name.c_str();
  unsigned int i = 0;
  for (unsigned int j = 1; j < name.size(); ++j) {
    if (c[j] == '[' || c[j] == ']' || c[j] == ',') {
      cut.push_back(name.substr(i,j-i));
      cut.push_back(name.substr(j,1));
      i = j+1;
    }
  }
  return cut;
}


static string makeMixtureName(MixtureNode const *mnode,
			      SymTab const &symtab)
{
  /* 
     Making a name for a mixture node is a complex business. We rely
     on the fact that, currently, mixture nodes can only be created by
     nested indexing. For example "foo[1,X,4]" where X can take values
     1,2,3, will create a mixture node with X as an index and
     parameters "foo[1,1,4]", "foo[1,2,4]", "foo[1,3,4]".

     The BUGS language name of a mixture node can be reconstructed
     by finding the parts of the parameter names that match, e.g.
     "foo[1," ... ",4]"
     and filling in the gaps with the names of the indices, e.g.
     "foo[1," "X" ",4]"
  */

  //The indices are the first elements of the vector of parents 
  vector<Node*> index = mnode->parents();
  index.resize(mnode->index_size());

  //Create a set of parents, excluding index nodes
  set<Node*> parents;
  for (unsigned int i = mnode->index_size(); i < mnode->parents().size(); ++i)
    {
      parents.insert(mnode->parents()[i]);
    }
  
  if (parents.empty())
    return("<unknown mixture node>");

  //Create a vector of parameter names
  vector<string> param_names;
  for (set<Node*>::const_iterator i = parents.begin(); i != parents.end(); ++i)
    {
      param_names.push_back((*i)->name(symtab));
    }

  //Cut up the name of the first parameter
  string const &name0 = param_names[0];
  vector<string> cut_name0 = cutBUGSSubsetName(name0);

  //Now compare the parameter names with the first one, and mark the
  //parts that don't match
  unsigned int nparts = cut_name0.size();
  vector<bool> variable_part(nparts, false);

  for (unsigned int j = 1; j < param_names.size(); ++j) {
    vector<string> cut_name = cutBUGSSubsetName(param_names[j]);
    if (cut_name.size() != nparts) {
      /* 
	 Something went wrong here, but it's too late to throw an exception
	 because we normally only want to know a node name after catching
	 one.
      */
      return "";
    }
    for (unsigned int k = 0; k < nparts; ++k) {
      if (cut_name[k] != cut_name0[k]) {
	variable_part[k] = true;
      }
    }
  }
    
  //Check that number of variables parts is the same as the number
  //of indices of the mixture node
  unsigned int nvar = 0;
  for (unsigned int k = 0; k < nparts; ++k) {
    if (variable_part[k])
      ++nvar;
  }
  if (nvar != index.size()) {
    /* Something went wrong - see above */
    return "";
  }
   
  //Now we are ready
  string mnode_name;
  unsigned int l = 0;
  for (unsigned int k = 0; k < nparts; ++k) {
    if (variable_part[k]) {
      mnode_name.append(index[l++]->name(symtab));
    }
    else {
      mnode_name.append(cut_name0[k]);
    }
  }

  return mnode_name;
}


string SymTab::getName(Node const *node) const
{
  //Special rules for mixture nodes
  if (isMixture(node)) {
    return makeMixtureName(asMixture(node), *this);
  }

  map<string, NodeArray*>::const_iterator p;
  for (p = _varTable.begin(); p != _varTable.end(); ++p) {
    NodeArray *array = p->second;
    if (array->graph().contains(node)) {
      if (array->getRange(node) == array->range()) {
         return p->first;
      }
      else {
         return p->first + print(array->getRange(node));
      }
    }
  }

  return("");
}

unsigned int SymTab::nchain() const
{
  return _nchain;
}
