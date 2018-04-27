#ifndef COMPILER_H_
#define COMPILER_H_

#include <compiler/LogicalFactory.h>
#include <compiler/MixtureFactory.h>
#include <compiler/CounterTab.h>
#include <compiler/ObsFuncTab.h>
#include <distribution/DistTab.h>
#include <function/FuncTab.h>
#include <model/BUGSModel.h>
#include <graph/Graph.h>

#include <map>
#include <string>
#include <utility>
#include <list>
#include <utility>
#include <set>

namespace jags {

class ParseTree;
class SymTab;
class FuncTab;
class DistTab;
class NodeAlias;
class SimpleRange;
class Compiler;

typedef void (Compiler::*CompilerMemFn) (ParseTree const *);

enum CompilerMode {PERMISSIVE, STOP_ON_ERROR, COLLECT_UNRESOLVED,
		   CLEAN_UNRESOLVED};
/**
 * @short Creates a BUGSModel from a ParseTree
 */
class Compiler {
  BUGSModel &_model;
  CounterTab _countertab;
  std::map<std::string, SArray> const &_data_table;
  std::map<std::string, std::vector<bool> > _constant_mask;
  unsigned int _n_resolved, _n_relations;
  std::vector<bool> _is_resolved;
  CompilerMode _compiler_mode;
  int _index_expression;
  std::vector<Node*> _index_nodes;
  LogicalFactory _logicalfactory;
  MixtureFactory _mixfactory1;
  MixtureFactory _mixfactory2;
  std::map<std::string, std::vector<unsigned long> > _node_array_bounds;
  std::map<std::pair<std::string, Range>, std::set<unsigned long> > _umap;
  std::set<std::string> _lhs_vars;
  std::map<std::pair<std::vector<unsigned long>, std::vector<double> >,
      ConstantNode *> _cnode_map;
  
  Node *getArraySubset(ParseTree const *t);
  SimpleRange VariableSubsetRange(ParseTree const *var);
  Range CounterRange(ParseTree const *var);
  Node* VarGetNode(ParseTree const *var);
  Range getRange(ParseTree const *var,  SimpleRange const &default_range);

  void traverseTree(ParseTree const *relations, CompilerMemFn fun,
		    bool resetcounter=true, bool reverse=false);
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
  Node * constFromTable(ParseTree const *p);
  void addDevianceNode();
  Node *getConstant(double value, unsigned int nchain, bool observed);
  Node *getConstant(std::vector<unsigned long> const &dim, 
		    std::vector<double> const &value,
		    unsigned int nchain, bool observed);
  void getLHSVars(ParseTree const *rel);
public:
  bool indexExpression(ParseTree const *t, std::vector<unsigned long> &value);
  BUGSModel &model() const;
  Node * getParameter(ParseTree const *t);
  /**
   * @param model Model to be created by the compiler.
   *
   * @param data_table Data table, mapping a variable name onto a
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
  /**
   * The table for observable functions used by the compiler to substitute
   * a logical node for a stochastic node when required
   */
  static ObsFuncTab &obsFuncTab();
  MixtureFactory &mixtureFactory1();
  MixtureFactory &mixtureFactory2();
};

} /* namespace jags */

#endif /* COMPILER_H_ */

