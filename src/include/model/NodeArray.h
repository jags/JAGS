#ifndef NODE_ARRAY_H_
#define NODE_ARRAY_H_

#include <graph/Graph.h>
#include <sarray/Range.h>

#include <string>
#include <map>

class SArray;

/**
 * A NodeArray is a container class for nodes with a dimension
 * attribute.  The array can be tiled with nodes using the insert()
 * function.  Inserted nodes can be retrieved using the find()
 * function.  Arbitrary subsets of the NodeArray can be returned with
 * the getSubSet() function.
 */
class NodeArray {
  std::string _name;
  Range _range;
  unsigned int _nchain;
  Graph _graph;
  Node **_node_pointers;
  long *_offsets;
  std::map<Range, Node *> _generated_nodes;

  /* Forbid copying */
  NodeArray(NodeArray const &orig);
  NodeArray &operator=(NodeArray const &rhs);

  bool findActiveIndices(Index &ind, unsigned int k, Index const &lower, Index const &dim) const;
public:
  /**
   * Constructor. Creates a NodeArray with the given name and dimension
   */
  NodeArray(std::string const &name, Index const &dim, unsigned int nchain);
  ~NodeArray();
  /**
   * Inserts a node into the subset of the NodeArray defined by range.
   * The dimension of the node must conform with the given range.
   * The given range must not overlap any previously inserted node.
   *
   * The node is added to an internal graph. 
   *
   * @exception runtime_error
   */
  void insert(Node *node, Range const &range);
  /**
   * Determines whether the given range is empty of inserted nodes,
   * and hence whether it can be used as an argument to insert
   */
  bool isEmpty(Range const &range) const;
  /**
   * Returns a node corresponding to the given range.  The range must
   * have been used in a previous call to insert, otherwise a NULL
   * pointer is returned.
   */
  Node* find(Range const &range) const;
  /**
   * Returns an arbitrary subset of the NodeArray. If the range
   * corresponds to a previously inserted node, this will be
   * returned. Otherwise, an aggregate node will be created, and it
   * will be added to the graph of the NodeArray.  If the range is not
   * completely covered by inserted nodes, a NULL pointer will be
   * returned.
   */
  Node* getSubset(Range const &range);
  /**
   * Sets the values of the nodes in the array. 
   *
   * @param value SArray containing values to be used. The value
   * vector of the SArray may contain missing values.  If so, then the
   * part of the value vector corresponding to each inserted node must be
   * either completely missing (in which no action is taken) or
   * contain no missing values.
   *
   * @param value array with dimensions matching the NodeArray containing
   * values.
   * @param chain Index number of chain to which to apply values
   * @param Set of nodes whose values have been altered.
   */
  void setValue(SArray const &value, unsigned int chain,
		std::set<Node*> &setnodes);
  /**
   * Gets the values of selected nodes that have been inserted into
   * the array.
   *
   * @param value SArray to which values should be written.
   *
   * @param chain Index number of chain to read.
   *
   * @param condition.  Boolean function that returns true for nodes
   * whose values are to be read.
   */
  void getValue(SArray &value, unsigned int chain,
		bool (*condition)(Node const *)) const;
  /**
   * Set data. All chains are set to the same value. If a value is
   * given for an index with no node, then a new ConstantNode is created
   * with that value.
   */
  void setData(SArray const &value);
  /**
   * Returns the name of the node array
   */
  std::string const &name() const;
  /**
   * Returns the range of indices covered by the NodeArray
   */
  Range const &range() const;
  /**
   * Returns the graph associated with the NodeArray.  This graph
   * contains all the nodes that were either used in a call to
   * insert, or created during a call to getSubset.
   */
  Graph const &graph() const;
  /**
   * Returns the range corresponding to the given node, if it
   * belongs to the graph associated with the NodeArray. If it is
   * not in the graph, a NULL Range is returned.
   */
  Range getRange(Node const *node) const;
  /**
   * Returns the number of chains of the nodes stored in the array
   */
  unsigned int nchain() const;
};


#endif /* NODE_ARRAY_H */
