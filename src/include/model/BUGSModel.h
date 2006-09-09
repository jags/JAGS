#ifndef BUGS_MODEL_H_
#define BUGS_MODEL_H_

#include <vector>
#include <map>
#include <utility>
#include <string>

#include <model/Model.h>
#include <model/SymTab.h>

/**
 * For a BUGS model we need to store extra information with the monitor
 * This is only required to get output in CODA format, which is a bit
 * of a nuisance.
 */

/*
struct BUGSMonitor
{
  NodeArray const *array;
  Range *range;
  TraceMonitor const *monitor;
};
*/

/**
 * @short BUGS model 
 *
 * A BUGS model is a subclass of Model that contains a symbol table
 * which is used to store certain nodes in arrays.  The array format
 * gives a convenient way of looking up nodes by name.
 */
class BUGSModel : public Model
{
  SymTab _symtab;
  std::map<Node const *, std::pair<std::string, Range> > _node_map;
 public:
  BUGSModel(unsigned int nchain);
  
  SymTab &symtab();
  void addDevianceNode();
  Node *getNode(std::string const &name, Range const &range,
		std::string &message);
  bool setMonitor(std::string const &name, Range const &range, unsigned int thin,
                  std::string &message);
  /*
  void clearMonitor(std::string const &name, Range const &range);
  TraceMonitor const *getMonitor(string const &name, Range const &range);
  */
  void coda(std::vector<Node const*> const &nodes, std::ofstream &index,
	    std::vector<std::ofstream*> &output);
  void setParameters(std::map<std::string, SArray> const &param_table,
		     unsigned int chain);
};

#endif /* BUGS_MODEL_H_ */
