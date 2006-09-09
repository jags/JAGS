#ifndef SYMTAB_H_
#define SYMTAB_H_
 
#include <graph/NodeNameTab.h>
#include <model/NodeArray.h>
#include <sarray/SArray.h>

#include <string>
#include <vector>
#include <map>

/**
 * @short Symbol table
 *
 * The SymTab class stores the names of variables used in the BUGS
 * language representation of the model.
 *
 * @see NodeArray
 */
class SymTab : public NodeNameTab
{
  unsigned int _nchain;
  std::map<std::string, NodeArray*> _varTable;
  std::map<Node const*, std::string> _names;
public:
  /**
   * Constructs an empty symbol table
   */
  SymTab(unsigned int nchain);
  /**
   * Adds an array variable to the symbol table. This creates a
   * NodeArray object of the given dimension and associates it
   * with the name, so it can be retrieved with a call to getVariable.
   * If no dimension is given, the variable is assumed to be scalar.
   */
  void addVariable(std::string const &name, Index const &dim = Index(1));
  /**
   * Returns a pointer to the  NodeArray associated with the given
   * name, or a NULL pointer if there is no such NodeArray.
   */
  NodeArray *getVariable(std::string const &name) const;
  /**
   * Inserts a node into the symbol table with the given name
   * (which must correspond to a previously added variable)
   * and range (which must be a valid sub-range for the variable).
   */
  void insertNode(Node *node, std::string const &name, Range const &range);
  /**
   * Adds all the nodes contained in all the NodeArray objects
   * in the SymTab to the given vector
   */
  void getNodes(std::vector<Node*> &nodes);
  /**
   *
   */
  void writeData(std::map<std::string, SArray> const &data_table);
  /**
   * Write values from the data table to the NodeArrays in SymTab with
   * the same name. Non-variable descendants of the nodes whose values
   * are set are updated.
   *
   * @param data_table Data table from which results will be read
   *
   * @param observed Logical flag. If true, then existing nodes will
   * have their value fixed after being written, and if there is no
   * existing node corresponding to a value in the data table, a
   * constant node will be created.
   */
  void writeValues(std::map<std::string, SArray> const &data_table,
		   unsigned int chain, std::set<Node*> &setnodes);
  /**
   * Reads the current value of selected nodes in the symbol table and
   * writes the result to the data table. The selection is based on
   * a given boolean function that returns true for the selected nodes.
   *
   * @param data_table Data table to which results will be written.
   * New entries will be created for the selected nodes.  However, a
   * new entry is not created if, in the symbol table, all nodes
   * corresponding to the selection are missing. Existing entries in
   * the data table will be overwritten.
   *
   * @param condition Function that returns true if the values of a
   * Node are to be read, and false otherwise.
   */
  void readValues(std::map<std::string, SArray> &data_table, 
                  unsigned int chain,
		  bool (*condition)(Node const *)) const;
  /**
   * Returns the number of variables in the symbol table
   */
  unsigned int size() const;
  /**
   * Deletes all the variables in the symbol table
   */
  void clear();
  /**
   * Gets the BUGS language name of the node if it belongs to
   * any of the NodeArrays in the symbol table. Special rules for nested
   * indexing also allow the names of Mixture Nodes to be calculated.
   * If the node name is not found, an empty string is returned
   */
  std::string getName(Node const *node) const;
  /**
   * Returns the number of chains that the symbol table can hold
   */
  unsigned int nchain() const;
};

#endif /* SYMTAB_H_ */
