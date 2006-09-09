#include <config.h>
#include <graph/Graph.h>
#include <graph/GraphMarks.h>
#include <graph/Node.h>

#include <stdexcept>
#include <vector>
#include <set>

using std::vector;
using std::set;
using std::invalid_argument;
using std::logic_error;

Graph::Graph()
  : _nodes()
{
}


Graph::~Graph()
{
  for (set<Node*>::iterator p = _nodes.begin(); p != _nodes.end(); p++) {
    (*p)->unref();
  }
}

void Graph::add(Node *node)
{
  if (!node) {
    throw invalid_argument("Attempt to add null node to graph");
  }
  if (!this->contains(node)) {
    node->ref();
    _nodes.insert(node);
  }
}

void Graph::remove(Node *node)
{
  if (this->contains(node)) {
    _nodes.erase(node);
    node->unref();
  }
}

bool Graph::contains(Node const *node) const
{
  return  _nodes.find(const_cast<Node*>(node)) != _nodes.end();
}

unsigned int Graph::size() const
{
  return _nodes.size();
}

bool Graph::isClosed()
{
  /* Determine whether any nodes in the graph have children or
     parents outside the graph */

  for (set<Node*>::iterator i = _nodes.begin(); 
       i != _nodes.end(); i++) {
    
    /* Check parents */
    vector<Node*> const &parents = (*i)->parents();
    for (vector<Node*>::const_iterator j = parents.begin(); 
	 j != parents.end(); j++) 
      {
	if (!this->contains(*j)) {
	  return false;
	}
      }

    /* Check children */
    set<Node*> const &children = (*i)->children();
    for (set<Node*>::iterator k = children.begin(); 
	 k != children.end(); k++) 
      {
	if (!this->contains(*k)) {
	  return false;
	}
      }
  }
  return true;
}

bool Graph::isConnected()
{
  GraphMarks marks(*this);
 
  /* Start by taking an arbitrary node */
  set<Node*>::const_iterator i = _nodes.begin();
  Node const *anode = *i;

  /* Mark the node, its parents and ancestors */
  marks.mark(anode, 1);
  marks.markAncestors(anode, 1);
  marks.markDescendants(anode, 1);
  
  /* If the graph is connected, all the other nodes will be marked */
  for (++i; i != _nodes.end(); ++i) {
    if (marks.mark(*i) == 0) {
      return false;
    }
  }
  return true;
}


bool Graph::isTree()
{
  /* 
     Determine whether a graph is a tree
     All nodes in the graph have exactly one parent in the graph,
     except for the root node which has none. The root node is unique
  */
  
  bool foundroot = false;

  for (set<Node*>::iterator i = _nodes.begin(); 
       i != _nodes.end(); i++) {
    
    int nparents = 0;
    vector<Node*> const &parents = (*i)->parents();
    for (vector<Node*>::const_iterator p = parents.begin(); 
	 p != parents.end(); ++i) {
      
      if (contains(*p))
	nparents++;
    }
      
    switch (nparents) {
    case 0:
      if (foundroot)
	return false;
      else 
	foundroot = true;
      break;
    case 1:
      break;
    default:
      return false;
      break;
    }
  }
  return true;
}

/*
static bool findDescendant(Node const *node, Node const *target, 
			   Graph const *graph, GraphMarks &marks)
{
  //Determine whether "target" is a descendent of "node" within the graph,
  //   excluding nodes that are already marked.

  if (marks.mark(node) == 1)
    return false;

  for (set<Node*>::const_iterator i = node->children().begin(); 
       i != node->children().end(); i++) 
    {
      if (graph->contains(*i)) {
	if (*i == target)
	  return true;
	if (findDescendant(*i, target, graph, marks))
	  return true;
      }
    }

  marks.mark(node, 1);
  return false;
}

bool Graph::hasCycle()
{
  GraphMarks marks(*this);
  for (set<Node*>::const_iterator p = _nodes.begin(); p != _nodes.end(); ++p) 
    {
      if (findDescendant(*p, *p, this, marks))
	return true;
    }
  return false;
}
*/

bool Graph::hasCycle()
{
  /* 
     Recursively grows a an acyclic graph G, consisting of marked
     nodes. If G grows to the whole graph then it is acyclic.
     
     We start with nodes that have no children. On each iteration
     we add nodes whose children are all in G. Adding such nodes
     cannot create a cycle.
  */

  GraphMarks marks(*this);
  while (true) {
    long marked = 0;
    long unmarked = 0;
    for (set<Node*>::const_iterator p = _nodes.begin(); p != _nodes.end();
	 ++p)
      {
	if (marks.mark(*p) == 0) {
	  bool can_mark = true;
	  set<Node*> const &children = (*p)->children();
	  for(set<Node*>::const_iterator c = children.begin();
	      c != children.end(); ++c) 
	    {
	      if (marks.mark(*c) == 0) {
		can_mark = false;
		break;
	      }
	    }
	  if (can_mark) {
	    marks.mark(*p,1);
	    ++marked;
	  }
	  else {
	    ++unmarked;
	  }
	}
      }
    if (marked == 0) {
      return (unmarked > 0);
    }
  }
}

set<Node*> const &Graph::nodes() const
{
  return _nodes;
}

void Graph::getNodes(vector<Node*> &nodes) const
{
  for (set<Node*>::iterator p = _nodes.begin();
       p != _nodes.end(); p++)
    {
      nodes.push_back(*p);
    }
}

void Graph::getSortedNodes(vector<Node*> &sortednodes)
{
  /* 
     Return a vector of nodes whose ordering follows the partial
     ordering of the graph.  If a is after b then there is never a
     path from a to b.

     Preconditions: Graph is acyclic
                    Sortednodes is empty

     The method will return false if any of these conditions is violated.
  */

  if (!sortednodes.empty()) {
    throw logic_error("vector not empty in Graph::getSortedNodes");
  }

  set<Node*>::const_iterator i;
  GraphMarks marks(*this);
  /* Start by marking the founder nodes */
  unsigned int n = _nodes.size();
  for (i = _nodes.begin(); i != _nodes.end(); i++) {
    Node *node = (*i);
    vector<Node*> const &parents = node->parents();
    bool isFounder = true;
    for (vector<Node*>::const_iterator j = parents.begin(); 
	 j != parents.end(); j++) 
      {
	if (this->contains(*j)) {
	  isFounder = false;
	  break;
	}
      }
    if (isFounder) {
      marks.mark(node, 1);
      sortednodes.push_back(node);
      n--;
    }
  }
  
  /* 
     Mark the rest of the nodes by successive passes through
     the graph until all nodes are marked
  */
  while (n > 0) {
    bool loopcheck = true;
    for (i = _nodes.begin(); i != _nodes.end(); i++) {
      Node *node = *i;
      if (marks.mark(node) == 0) {

	/* We can only mark a node if its parents in the graph are marked */
	vector<Node*> const &parents = node->parents();
	bool can_mark = true;
	for (vector<Node*>::const_iterator j = parents.begin(); 
	     j != parents.end(); j++) 
	  {
	    Node *parent = *j;
	    if (this->contains(parent) && marks.mark(parent) == 0) {
	      /* Unmarked parent */
	      can_mark = false;
	      break;
	    }
	  }
	   
	if (can_mark) {
	  marks.mark(node, 1);
	  sortednodes.push_back(node);
	  n--;
	  loopcheck = false;
	}
      }
    }
    if (loopcheck) {
      throw logic_error("Failure in Graph::getSortedNodes. Directed cycle in graph");
    }
  }
}
