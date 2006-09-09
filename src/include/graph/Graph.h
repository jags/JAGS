#ifndef GRAPH_H_
#define GRAPH_H_

#include <set>
#include <vector>
class Node;

/**
 * A graph is a container class for (pointers to) Nodes. A Node may
 * belong to several Graphs. Further, if Node N is in graph G, then
 * there is no requirement that the parents or children of N lie in G.
 *
 * @short Container class for nodes
 */
class Graph {
  std::set<Node*> _nodes;

  /* forbid copying */
  Graph(Graph const &orig);
  Graph &operator=(Graph const &rhs);
public:
  /**
   * Constructs an empty graph.
   */
  Graph();
  /**
   * Destructor.  The reference count of all nodes in the graph
   * is decremented.
   */
  ~Graph();
  /**
   * Adds node to graph.  The reference count of the node is
   * incremented.  If node is already in the graph, no action is
   * taken.
   */
  void add(Node *node);
  /**
   * Removes node from graph. The reference count of the node is
   * decremented. If node is not a member, no action is taken
   */
  void remove(Node *node);
  /**
   * Checks to see whether the node is contained in the Graph.
   */
  bool contains(Node const *node) const;
  /**
   * The number of nodes in the graph.
   */
  unsigned int size() const;
  /**
   * Checks if the graph is connected.
   */
  bool isConnected();
  /**
   * Checks if the parents and children of every node in the
   * graph are also contained in the graph.
   */
  bool isClosed();
  /**
   * Checks if there is any path in the graph leading from a
   * node to itself.
   */
  bool hasCycle();
  /**
   * Checks if every node in the graph has exactly one parent in
   * the graph except for a unique root node which has none.
   */
  bool isTree();
  /**
   * The set of nodes contained in the graph
   */
  std::set<Node*> const &nodes() const;
  /**
   * Adds all nodes in the graph to the given vector
   */
  void getNodes(std::vector<Node*> &nodes) const;
  /**
   * Adds all nodes in the graph to the given vector with partial
   * ordering, so that if A is an ancestor of B, then B never appears
   * before A in the vector (Note that if there is a path from A to B
   * outside of the graph, then this is ignored).
   *
   * The graph must be acyclic.
   *
   * @param sorted Empty vector of Node pointers.  On exit
   * this vector will contain the sorted nodes.
   */
  void getSortedNodes(std::vector<Node*> &sorted); 
};

#endif /* GRAPH_H_ */
