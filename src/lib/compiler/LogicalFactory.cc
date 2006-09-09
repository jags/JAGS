#include <config.h>
#include <compiler/LogicalFactory.h>
#include <graph/Node.h>
#include <graph/LogicalNode.h>
#include <function/Function.h>

#include <stdexcept>

using std::pair;
using std::map;
using std::vector;
using std::invalid_argument;

static bool equal(vector<Node*> const &v1, vector<Node*> const &v2)
{
  unsigned long n1 = v1.size();
  unsigned long n2 = v1.size();
  if (n1 != n2) {
    return false;
  }
  else {
    for (unsigned int i = 0; i < n1; ++i) {
      if (!equal(v1[i], v2[i])) {
	return false;
      }
    }
    return true;
  }
}

bool lt(LogicalPair const &arg1, LogicalPair const &arg2)
{
  // Order first by function
  if (arg1.first < arg2.first) {
    return true;
  }
  else if (arg1.first > arg2.first) {
    return false;
  }
  else {
    //Same function. Check parameters
    //Order by parameter length
    unsigned long n1 = arg1.second.size();
    unsigned long n2 = arg2.second.size();
    if (n1 < n2) {
      return true;
    }
    else if (n1 > n2) {
      return false;
    }
    else {
      //Same parameter length. 
      //Check for equality
      if (equal(arg1.second, arg2.second)) {
	return false;
      }
      else {
	//Order by address in memory
	for (unsigned int i = 0; i < n1; ++i) {
	  if (arg1.second[i] < arg2.second[i]) {
	    return true;
	  }
	  else if (arg1.second[i] > arg2.second[i]) {
	    return false;
	  }
	}
	return false;
      }
    }
  }
}

LogicalNode*
LogicalFactory::getLogicalNode(Function const *func, 
			       vector<Node*> const &parents)
{
  if (func == 0) {
    throw invalid_argument("NULL function passed to getLogicalNode");
  }

  LogicalPair lpair(func, parents);
  map<LogicalPair, LogicalNode*>::iterator i = _logicalmap.find(lpair);

  if (i != _logicalmap.end()) {
    return i->second;
  }
  else {
    // Create a new logical node
    LogicalNode *lnode = new LogicalNode(func, parents);
    _logicalmap[lpair] = lnode;
    _graph.add(lnode);
    return lnode;
  }
}
