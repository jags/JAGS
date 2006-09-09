#ifndef COUNTERTAB_H_
#define COUNTERTAB_H_

#include <model/NodeArray.h>
#include <compiler/Counter.h>

#include <string>
#include <vector>
#include <map>

/**
 * @short Counter table
 *
 * The CounterTab class stores counters in the BUGS language
 * representation of the model. Counters are only required during
 * compilation of the model.
 *
 * @see NodeArray
 * @see Counter
 */
class CounterTab {
  std::vector<std::pair<std::string, Counter*> > _table;
public:
  CounterTab();
  ~CounterTab();
  /**
   * Pushes a counter onto the internal counter stack and returns
   * a pointer to it.
   * @param name Name of the counter
   * @param range Range over which the counter will vary
   */
  Counter * pushCounter(std::string const &name, Range const &range);
  /**
   * Pops the top counter off the counter stack, freeing the associated
   * memory in the process.
   */
  void popCounter();
  /**
   * Returns a pointer associated with the given name, or a NULL
   * pointer if there is no such Counter
   */
  Counter *getCounter(std::string const &name) const;
};

#endif /* COUNTERTAB_H_ */

