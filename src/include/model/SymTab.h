#ifndef SYMTAB_H_
#define SYMTAB_H_
 
#include <model/NodeArray.h>
#include <sarray/SArray.h>

#include <string>
#include <vector>
#include <map>

namespace jags {

/**
 * @short Associates a NodeArray object with a BUGS-language name
 *
 * The SymTab class stores the names of variables used in the BUGS
 * language representation of the model.
 *
 * @see NodeArray
 */
class SymTab  
{
  Model *_model;
  std::map<std::string, NodeArray*> _varTable;
  std::map<Node const*, std::string> _names;
public:
  /**
   * Constructs an empty symbol table
   *
   * @param model Model to which newly allocated nodes are added.
   *
   */
  SymTab(Model *model);
  ~SymTab();
  /**
   * Adds an array variable to the symbol table. This creates a
   * NodeArray object of the given dimension and associates it
   * with the name, so it can be retrieved with a call to getVariable.
   * If no dimension is given, the variable is assumed to be scalar.
   */
  void addVariable(std::string const &name, std::vector<unsigned long> const &dim);
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
   * Creates constant nodes in all the arrays of the symbol table
   * taking values from the given data table.
   *
   * @param data_table Data table from which values will be read.
   *
   * @see NodeArray#setData
   */
  void writeData(std::map<std::string, SArray> const &data_table);
  /**
   * Write values from the data table to the arrays in the symbol
   * table with the same name. Unlike the writeData function, the
   * Nodes are not permanently set to the supplied values, and values
   * are only written to the given chain.
   *
   * @param data_table Data table from which values will be read
   *
   * @param chain Index number of chain to which values are written.
   *
   * @see NodeArray#setValue
   */
  void writeValues(std::map<std::string, SArray> const &data_table,
		   unsigned int chain);
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
  unsigned long size() const;
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
};

} /* namespace jags */

#endif /* SYMTAB_H_ */
