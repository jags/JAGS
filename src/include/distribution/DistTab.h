#ifndef DIST_TAB_H_
#define DIST_TAB_H_

#include <map>
#include <string>

class Distribution;

/** 
 * A DistTab is a table which can be used to look up distributions by
 * name. A DistTab does no memory management: the module responsible
 * for creating any Distribution inserted into the table is also
 * responsible for deleting it.
 *
 * @short Distribution table 
 * 
 * @see Distribution
 */
class DistTab {
  std::map<const std::string, const Distribution*> _table;
public:
  /**
   * Distributions in the table should be removed before the DistTab
   * is deleted.  The destructor therefore prints a warning to the
   * standard error if the DistTab contains any Distributions when it
   * is called.
   */
  ~DistTab();
  /**
   * Inserts a distribution into the table.  
   * 
   * @return logical value indicating success.  The insertion will
   * fail if a distribution with the same name is already in the
   * DistTab.  
   *
   * @param dist Pointer to the distribution to insert, which should
   * be dynamically allocated.
   */
  bool insert(Distribution const *dist);
  /** 
   * Removes a distribution by name and frees the associated memory.
   * 
   * With a combination of erase and insert calls, you can replace any
   * existing distribution with another implementation.
   *
   * @return a logical value indicating success.  The removal will
   * fail if there is no distribution with the given name in the
   * DistTab.
   */
  bool erase(std::string const &name);
  /**
   * Finds a distribution by name
   *
   * @return a pointer to the distribution or 0 if it was not found.
   */
  Distribution const *find(std::string const &name) const;
  /**
   * Removes a distribution, if it is contained in the table: does
   * nothing if the distribution is not in the table.
   *
   * In practice, this function can only be called by the module that
   * created the distribution, since it requires a non-constant
   * pointer
   */
  void erase(Distribution *dist);
};

#endif /* DIST_TAB_H_ */
