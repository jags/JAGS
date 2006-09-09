#ifndef FUNC_SORT_H_
#define FUNC_SORT_H_

#include <function/Function.h>

/**
 * @short Sort function
 * Sorts the elements of a vector in ascending order
 * <pre>
 * y <- sort(x[])
 * </pre>
 */
class Sort : public Function
{
public:
  Sort ();
  void evaluate (SArray & value, std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  Index dim (std::vector <Index> const &dims) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
};

#endif /* FUNC_SORT_H_ */
