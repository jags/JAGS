#ifndef GRAPH_MARKS_H_
#define GRAPH_MARKS_H_

#include <map>
class Node;
class Graph;

/**
 * A GraphMarks object associates each node in the graph with an
 * integer mark.  This is useful for several algorithms.  Rather than
 * modifying the original graph, the marks are kept in a separate
 * table.
 *
 * @short Mark a graph.
 * @see Graph
 */
class GraphMarks {
  std::map<Node const*,int> _marks;
 public:
  /** Constructor. Each node in the graph initially has mark zero */
  GraphMarks(Graph const &graph);
  ~GraphMarks();

  /** Sets the mark of node to m. The graph must contain node */
  void mark(Node const *node, int m);
  /** Gives the mark of node */
  int mark(Node const *node);
  /** Sets the mark all nodes in the graph to m */
  void markAll(int m);
  /** Marks the parents of node that are in the graph.
      The graph must contain node */
  void markParents(Node const *node, int m);
  /** Marks the children of node that are in the graph. The
      graph must contain node */
  void markChildren(Node const *node, int m);
  /** Marks the ancestors of the node in the graph, i.e. every
      node N for which there is a directed path from N to node
      in the graph. The graph must contain node
  */
  void markAncestors(Node const *node, int m);
  /** Marks the descendants of the node in the graph, i.e. every
      node N for which there is a directed path from node N
      in the graph. The graph must contain node.
  */
  void markDescendants(Node const *node, int m);

};

#endif /* GRAPH_MARKS_H_ */
