#include <config.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <graph/Node.h>
#include <graph/NodeError.h>

#include <vector>
#include <set>

using std::map;
using std::vector;
using std::set;
using std::pair;

GraphMarks::GraphMarks(Graph const &graph)
{
  set<Node*> const &nodes = graph.nodes();
  for (set<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i) {
    _marks[*i] = 0;
  }
}

GraphMarks::~GraphMarks()
{}

void GraphMarks::mark(Node const *node, int m)
{
  map<Node const*, int>::iterator i = _marks.find(node);
  if (i == _marks.end()) {
    throw NodeError(node, "Attempt to set mark of node not in Graph");
  }
  else {
    i->second = m;
  }
}

int GraphMarks::mark(Node const *node)
{
  map<Node const*, int>::iterator i = _marks.find(node);
  if (i == _marks.end()) {
    throw NodeError(node, "Attempt to get mark of node not in Graph");
  }
  return i->second;
}

void GraphMarks::markAll(int m)
{
  map<Node const*, int>::iterator i = _marks.begin();
  for (;i != _marks.end(); ++i) {
    _marks[i->first] = m;
  }
}

void GraphMarks::markParents(Node const *node, int m)
{

  if (_marks.find(node) == _marks.end()) {
    throw NodeError(node, "Can't mark parents of node: not in Graph");
  }
  else {
    vector<Node*> const &parents = node->parents();
    for (vector<Node*>::const_iterator p = parents.begin(); p != parents.end(); 
	 ++p) 
      {
	map<Node const*, int>::iterator i = _marks.find(*p);    
	if (i != _marks.end()) {
	  i->second = m;
	}
      }
  }
}
 
void GraphMarks::markChildren(Node const *node, int m)
{

  if (_marks.find(node) == _marks.end()) {
    throw NodeError(node, "Can't mark children of node: not in Graph");
  }
  else {
    set<Node*> const &children = node->children();
    for (set<Node*>::const_iterator p = children.begin(); p != children.end(); 
	 ++p) 
      {
	map<Node const*, int>::iterator i = _marks.find(*p);    
	if (i != _marks.end()) {
	  i->second = m;
	}
      }
  }
}

void GraphMarks::markDescendants(Node const *node, int m)
{
  /* 
     We need a temporary set of marks in order to keep track of which
     nodes in the graph have already been visited. This prevents an
     infinite loop in the case where the graph has a cycle.

     FIXME: Yes but memory usage is a disaster.

  GraphMarks tempmarks(*this);
  tempmarks.markAll(0);

  for (set<Node*>::const_iterator i = node->children().begin();
       i != node->children().end(); ++i) {
    if (_marks.find(*i) != _marks.end() && tempmarks.mark(*i) == 0) {
      mark(*i, m);
      tempmarks.mark(*i, 1);
      markDescendants(*i, m);
    }
  }
  */
  for (set<Node*>::const_iterator i = node->children().begin();
       i != node->children().end(); ++i) {
    if (_marks.find(*i) != _marks.end()) {
      mark(*i, m);
      markDescendants(*i, m);
    }
  }
}

void GraphMarks::markAncestors(Node const *node, int m)
{
  /* See markDescendants */
  /*
  GraphMarks tempmarks(*this);
  tempmarks.markAll(0);

  for (set<Node*>::const_iterator i = node->parents().begin();
       i != node->parents().end(); ++i) {
    if (_marks.find(*i) != _marks.end() && tempmarks.mark(*i) == 0) {
      mark(*i, m);
      tempmarks.mark(*i, 1);
      markAncestors(*i, m);
    }
  }
  */
  for (vector<Node*>::const_iterator i = node->parents().begin();
       i != node->parents().end(); ++i) {
    if (_marks.find(*i) != _marks.end()) {
      mark(*i, m);
      markAncestors(*i, m);
    }
  }
}
