#ifndef FUNC_RANK_H_
#define FUNC_RANK_H_

#include <function/Function.h>

/**
 * @short Sort function
 * Calculates the ranks of the elements of a vector in ascending order
 * <pre>
 * y <- sort(x[])
 * </pre>
 */
class Rank : public Function
{
public:
  Rank ();
  void evaluate (SArray & value, std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  Index dim (std::vector <Index> const &dims) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
};

#endif /* FUNC_RANK_H_ */
