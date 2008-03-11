#ifndef COMPILER_H_
#define COMPILER_H_

#include <compiler/LogicalFactory.h>
#include <compiler/ConstantFactory.h>
#include <compiler/MixtureFactory.h>
#include <compiler/CounterTab.h>
#include <distribution/DistTab.h>
#include <function/FuncTab.h>
#include <model/BUGSModel.h>

#include <map>
#include <string>
#include <utility>
#include <list>

class ParseTree;
class Graph;
class SymTab;
class FuncTab;
class DistTab;
class NodeAlias;

class Compiler;
typedef void (Compiler::*CompilerMemFn) (ParseTree const *);

/**
 * @short Creates a BUGSModel from a ParseTree
 */
class Compiler {
  BUGSModel &_model;
  CounterTab _countertab;
  std::map<std::string, SArray> const &_data_table;
  std::map<std::string, std::vector<bool> > _constant_mask;
  unsigned int _n_resolved, _n_relations;
  bool *_is_resolved;
  bool _strict_resolution;
  int _index_expression;
  Graph _index_graph;
  ConstantFactory _constantfactory;
  LogicalFactory _logicalfactory;
  MixtureFactory _mixfactory;
  std::map<std::string, std::vector<std::vector<int> > > _node_array_ranges;

  Node *getArraySubset(ParseTree const *t);
  Range VariableSubsetRange(ParseTree const *var);
  Range CounterRange(ParseTree const *var);
  Node* VarGetNode(ParseTree const *var);
  Range getRange(ParseTree const *var,  Range const &default_range);

  void traverseTree(ParseTree const *relations, CompilerMemFn fun,
		    bool resetcounter=true);
  void allocate(ParseTree const *rel);
  Node * allocateStochastic(ParseTree const *stoch_rel);
  Node * allocateLogical(ParseTree const *dtrm_rel);
  void setConstantMask(ParseTree const *rel);
  void writeConstantData(ParseTree const *rel);
  Node *getLength(ParseTree const *p, SymTab const &symtab);
  Node *getDim(ParseTree const *p, SymTab const &symtab);
  void getArrayDim(ParseTree const *p);
  bool getParameterVector(ParseTree const *t,
			  std::vector<Node const *> &parents);
  Node * getSubsetNode(ParseTree const *var);
  Node * constFromTable(ParseTree const *p);
  void addDevianceNode();
public:
  bool indexExpression(ParseTree const *t, int &value);
  BUGSModel &model() const;
  Node * getParameter(ParseTree const *t);
  /**
   * @param model Model to be created by the compiler.
   *
   * @param datatab Data table, mapping a variable name onto a
   * multi-dimensional array of values. This is required since some
   * constant expressions in the BUGS language may depend on data
   * values.
   */
  Compiler(BUGSModel &model, std::map<std::string, SArray> const &data_table);
  /**
   * Adds variables to the symbol table.
   *
   * @param pvariables vector of ParseTree pointers, each one corresponding
   * to a parsed variable declaration.
   */
  void declareVariables(std::vector<ParseTree*> const &pvariables);
  /**
   * Adds variables without an explicit declaration to the symbol
   * table.  Variables supplied in the data table are added, then
   * any variables that appear on the left hand side of a relation
   *
   * @param prelations ParseTree corresponding to a parsed model block
   */
  void undeclaredVariables(ParseTree const *prelations);
  /**
   * Traverses the ParseTree creating nodes.
   *
   * @param prelations ParseTree corresponding to a parsed model block
   */
  void writeRelations(ParseTree const *prelations);
  /**
   * The function table used by the compiler to look up functions by
   * name.  It is shared by all Compiler objects.
   *
   * @see Module
   */
  static FuncTab &funcTab();
  /**
   * The distribution table used by the compiler to look up
   * distributions by name.  It is shared by all Compiler objects.
   *
   * @see Module
   */
  static DistTab &distTab();
  MixtureFactory &mixtureFactory();
};


#endif /* COMPILER_H_ */

