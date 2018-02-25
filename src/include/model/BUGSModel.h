#ifndef BUGS_MODEL_H_
#define BUGS_MODEL_H_

#include <vector>
#include <map>
#include <utility>
#include <string>

#include <model/Model.h>
#include <model/SymTab.h>
#include <model/MonitorInfo.h>

namespace jags {

/**
 * @short Model with symbol table 
 *
 * A BUGS model is a subclass of Model that contains a symbol table
 * which is used to store certain nodes in arrays.  The array format
 * gives a convenient way of looking up nodes by name.
 */
class BUGSModel : public Model
{
    SymTab _symtab;
    //std::map<Node const*, std::pair<std::string, Range> > _node_map;
    std::list<MonitorInfo> _bugs_monitors;

    // Only to be used by observedStochasticNodes():
    std::vector<Node const *> _observed_stochastic_nodes;
    // Only to be used by observedStochasticNodeNames():
    std::vector<std::string> _observed_stochastic_node_names;

public:
    BUGSModel(unsigned int nchain);
    ~BUGSModel();
    /**
     * Returns the symbol table of the BUGSModel.
     */
    SymTab &symtab();
    /**
     * Writes out selected monitors in CODA format.
     *
     * @param nodes Vector of nodes to write out. Each node is
     * described by a pair consisting of a name and a range of
     * indices.  If a node is not being monitored, then it is ignored.
     *
     * @param prefix String giving prefix to be prepened to the output
     * file names.
     *
     * @param warn String that will contain any warning messages on
     * exit. It is cleared on entry.
     * 
     * @param type Name of the monitor type or "*" for all types
     *
     * @exception logic_error
     */
    void coda(std::vector<std::pair<std::string,Range> > const &nodes, 
	      std::string const &prefix, std::string &warn, std::string const &type);
    /**
     * Write out all monitors in CODA format
     */
    void coda(std::string const &prefix, std::string &warn, std::string const &type);
    /**
     * Sets the state of the RNG, and the values of the unobserved
     * stochastic nodes in the model, for a given chain.
     *
     * @param param_table STL map, in which each entry relates a
     * variable name to an SArray.  If the name is ".RNG.state" or
     * ".RNG.seed", then the SArray is used to set the state of the RNG.
     * Otherwise the SArray value is used to set the unobserved
     * stochastic nodes in the model. Each SArray must have the same
     * dimensions as the NodeArray in the symbol table with the
     * corresponding name. Elements of the SArray must be set to the
     * missing value, JAGS_NA, unless they correspond to an element of
     * an unobserved StochasticNode.
     *
     * @param chain Number of chain (starting from zero) for which
     * parameter values should be set.
     *
     * @see RNG#init RNG#setState
     * @exception runtime_error
     */
    void setParameters(std::map<std::string, SArray> const &param_table,
		       unsigned int chain);
    /**
     * Creates a new Monitor. The BUGSModel is responsible for the
     * memory management of any monitor created this way. It is not
     * possible to create two monitors with the same name, range and
     * type.
     *
     * @param name Name of the node array
     *
     * @param range Subset of indices of the node array defining hte
     * node to be monitored.
     * 
     * @param thin Thinning interval for monitor
     *
     * @param type Type of monitor to create
     *
     * @param msg User-friendly error message that may be given if no
     * monitor can be created.
     *
     * @return True if the monitor was created.  
     */
    bool setMonitor(std::string const &name, Range const &range,
		    unsigned int thin, std::string const &type,
		    std::string &msg);
    /**
     * Deletes a Monitor that has been previously created with a call
     * to setMonitor.
     *
     * @return True if the monitor was deleted.
     */
    bool deleteMonitor(std::string const &name, Range const &range,
		       std::string const &type);
    /**
     * Traverses the list of monitor factories requesting default
     * monitors of the given type. The function returns true after the
     * first monitor factory has added at least one node to the monitor
     * list. If none of the available monitor factories can create
     * default monitors of the given type, the return value is false.
     *
     * @see MonitorFactory#addDefaultMonitors
     */
    bool setDefaultMonitors(std::string const &type, unsigned int thin);
    /**
     * Removes all Monitors of the given type.
     */
    void clearMonitors(std::string const &type);
    /**
     * Writes the names of the samplers, and the corresponding 
     * sampled nodes vectors to the given vector.
     *
     * @param sampler_names vector that is modified during the call On
     * exit it will contain an element for each Sampler in the model.
     * Each element is a vector of strings: the first string is the
     * name of the sampler, and the remaining strings are the names of
     * the nodes sampled by that Sampler.
     */
    void samplerNames(std::vector<std::vector<std::string> > &sampler_names) 
	const;

    /**
     * Returns a vector of all observed stochastic nodes in the model
     */ 
    std::vector<Node const *> const &observedStochasticNodes();
    /**
     * Returns a vector of the names corresponding to the observed
     * stochastic nodes in the model
     */ 
    std::vector<std::string> const &observedStochasticNodeNames();
};

} /* namespace jags */

#endif /* BUGS_MODEL_H_ */
