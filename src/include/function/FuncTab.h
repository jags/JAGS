#ifndef FUNC_TAB_H_
#define FUNC_TAB_H_

#include <list>
#include <string>

class InverseLinkFunc;
class Function;

/**
 * @short Look-up table for Function objects
 *
 * Since all member functions of the JAGS Function class are constant, only
 * one instance of a Function object is required. The FuncTab class 
 * provides a means of storing Functions and looking them up by name.
 *
 * @see Function DistTab 
 */
class FuncTab
{
    std::list<Function const *> _func_list;
    std::list<Function const *> _masked_func_list;
    std::list<InverseLinkFunc const *> _link_list;
    std::list<InverseLinkFunc const *> _masked_link_list;
public:
  /**
   * Inserts an inverse link function into the table. This function
   * works the same way as the standard insert function (see below)
   * but additionally registers the function by its link name.
   */
  void insert (InverseLinkFunc const *func);
  /**
   * Inserts a function into the table. This function will mask any
   * function of the same name previously inserted into the table,
   * making them inaccessible until this function is erased.
   *
   * @param func Function to insert into the table.
   */
  void insert (Function const *func);
  /**
   * Finds a function by name.  If more than one function with the
   * same name has been inserted into the table, the most recently
   * inserted Function will be returned.
   *
   * @return a pointer to the function or a NULL pointer if it was not
   * found.
   */
  Function const *find (std::string const &name) const;
  /**
   * Finds the inverse of a link function by the link name
   *
   * @return a pointer to the inverse function or a NULL pointer if it
   * was not found.
   */
  Function const *findInverse (std::string const &name) const;
  /**
   * Finds an inverse link function by its name.  This is a more
   * restrictive version of the FuncTab#find function which returns
   * a pointer to an InverseLinkFunc object, or a NULL pointer if the
   * name does not correspond to an inverse link function.
   */
  InverseLinkFunc const *findLink (std::string const &name) const;
  /**
   * Removes an inverse link function from the table
   */
  void erase(InverseLinkFunc *func);
  /**
   * Removes a function from the table. This can only be called by the
   * module that owns the Function object, as it requires a
   * non-constant pointer.
   */
  void erase(Function *func);
};

#endif /* FUNC_TAB_H_ */
