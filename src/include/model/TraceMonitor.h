#ifndef TRACE_MONITOR_H_
#define TRACE_MONITOR_H_

class Node;

/**
 * A trace monitor stores every single sampled value of a given
 * node.
 */
class TraceMonitor {
  Node const *_node;                  
  long _start;                  
  long _end;                    
  int _thin;                         // thinning interval (=1/frequency)
  double *_values;                   // vector containing sampled values
  long _size;                        // current size of _values
  long _maxsize;                     // maximum size of _values (or -1)
  long _current;                     // current index value of vector
 public:
  /**
   * Constructor for a TraceMonitor of finite size
   * @param node Node to sample
   * @param start Iteration at which to start monitoring
   * @param end Iteration at which to end monitoring
   * @param thin Thinning interval for monitor.  The monitor only
   * stores iteration for which (iteration number)%thin == 0.
   */
  TraceMonitor(Node const *node, long start, long end, int thin);
  /**
   * Constructor for an open-ended TraceMonitor. The Monitor will
   * continue to sample values until it is deleted
   * @param node Node to sample
   * @param start Iteration at which to start monitoring
   * @param thin Thinning interval for monitor.  The monitor only
   * stores iterations for which (iteration number)%thin == 0.
   */
  TraceMonitor(Node const *node, long start, int thin);
  ~TraceMonitor();
  /**
   * Iteration number at which the node started monitoring.  This
   * may be different from the starting iteration requested in the
   * constructor if the thinning interval is > 1.
   */
  long start() const; 
  /**
   * For a finite monitor, returns the iteration number at which the
   * node will stop monitoring.  For an open ended monitor, returns
   * the last monitored node.
   */
  long end() const;
  /**
   * Thinning interval of monitor
   */
  long thin() const;
  /**
   * Sample size: the number of monitored iterations
   */
  long size() const;                       
  /**
   * The number of monitored iterations between two limits
   */
  long size(long start, long end) const;
  /**
   * Updates the monitor
   * @param iteration The current iteration number.
   */
  void update(long iteration, unsigned int chain);
  /**
   * Returns a pointer to the start of the vector of monitored values.
   * For an open-ended monitor, his pointer may be invalidated by a
   * subsequent call to update, since the value vector may be reallocated
   */
  double const *values() const; 
  /** 
   * Returns the monitored node
   */
  Node const *node() const;
};

#endif
