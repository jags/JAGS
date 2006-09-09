#ifndef FUNC_TAB_H_
#define FUNC_TAB_H_

#include <map>
#include <string>

class Function;
class InverseLinkFunc;

/** 
 * The FuncTab class provides a convenient way of storing Function
 * objects and referencing them by name.
 *
 * @short Function table
 * @see Function
 */
class FuncTab
{
    std::map < const std::string, const Function *> _table;
    std::map < const std::string, const InverseLinkFunc *> _link_table;
public:
  /**
   * Any functions inserted into the FuncTab should be removed before
   * it is destroyed. The destructor therefore prints a warning
   * message to the standard error if the table is non-empty when it
   * is called
   */
  ~FuncTab ();
  /**
   * Inserts an inverse link function into the table. This function
   * works the same way as the standard insert function (see below)
   * but additionally registers the function by its link name 
   * @see InverseLinkFunc
   */
  bool insert (InverseLinkFunc const *func);
  /**
   * Inserts a function into the table.
   *
   * @return logical value indicating success.  The insertion will
   * fail if a function with the same name is already in the FuncTab.
   *
   * @param func Pointer to the function to insert, which should
   * be dynamically allocated. The FuncTab is responsible for freeing
   * the memory allocated to func.
   */
  bool insert (Function const *func);
  /** 
   * Removes a function by name
   * 
   * With a combination of erase and insert calls, you can
   * replace any existing function with another implementation.
   *
   * @return a logical value indicating success.  The removal will
   * fail if there is no function with the given name in the FuncTab.
   */
  bool erase (std::string const &name);
  /**
   * Finds a function by name 
   *
   * @return a pointer to the function or 0 if it was not found.
   */
  Function const *find (std::string const &name) const;
  /**
   * Finds the inverse of a link function by the link name
   *
   * @return a pointer to the inverse function or 0 if it was not found.
   */
  Function const *findInverse (std::string const &name) const;
  /**
   * Removes a function from the table. This can only be called by the
   * module that owns the Function object, as it requires a
   * non-constant pointer.
   */
  void erase(Function *func);
};

#endif /* FUNC_TAB_H_ */
