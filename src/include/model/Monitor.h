#ifndef MONITOR_H_
#define MONITOR_H_

#include <sarray/SArray.h>

#include <vector>
#include <string>

class Node;

/**
 * @short Analyze sampled values 
 *
 * This is an abstract class for objects that analyze and/or store sampled
 * values from a given node. 
 */
class Monitor {
    std::string _type;
    Node const * _node;
    unsigned int _start;
    unsigned int _thin;
    unsigned int _niter;
public:
    Monitor(std::string const &type, Node const *nodes, unsigned int start,
            unsigned int thin);
    virtual ~Monitor();
    /**
     * Updates the monitor. If the iteration number coincides with
     * the thinning interval, then the doUpdate function is called.
     *
     * @param iteration The current iteration number.
     */
    void update(unsigned int iteration);
    /**
     * Returns the sampled Node
     */
    Node const *node() const;
    /**
     * @returns the iteration number at which the node started monitoring.  
     */
    unsigned int start() const; 
    /**
     * @returns The last monitored iteration
     */
    unsigned int end() const;
    /**
     * @returns the thinning interval of the monitor
     */
    unsigned int thin() const;
    /**
     * The number of monitored iterations
     */
    unsigned int niter() const;      
    /**
     * The type of monitor. Each subclass must define have a unique
     * type, which is common to all Monitors of that class. The type
     * is used by the user-interface to identify the subclass of Monitor.
     */
    std::string const &type() const;
    /**
     * The number of parallel chains of the Monitor.  This does not
     * have to coincide with the number of chains of the Model: a
     * Monitor may summarize data from multiple parallel chains in a
     * single vector.
     */
    virtual unsigned int nchain() const  = 0;
    /**
     * The dimension of the value vector for a single chain.
     */
    virtual std::vector<unsigned int> dim() const = 0;
    /**
     * The vector of monitored values for the given chain
     */
    virtual std::vector<double> const &value(unsigned int chain) const = 0;
    /**
     * Updates the Monitor with Node valus from the current iteration
     */
    virtual void doUpdate() = 0;
    /**
     * Reserves memory for future updates. Sufficient memory is
     * reserved for storage of future samples to avoid re-allocation
     * of memory for the next niter iterations.
     *
     * @param niter number of future iterations to reserve. 
     */
    virtual void reserve(unsigned int niter) = 0;
     /**
      * Dumps the monitored values to an SArray. 
      *
      * The SArray should have informative dimnames. In particular, the
      * dimnames "iteration" and "chain" should be used if there are
      * values for each iteration and each chain, respectively.
      */
     virtual SArray dump() const = 0;
};

#endif
