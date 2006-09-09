#ifndef FUNC_SUM_H_
#define FUNC_SUM_H_

#include <function/ScalarFunc.h>

/**
 * @short Sum of an array
 * Sum calculates the sum of the elements of an array
 * @see Mean
 * <pre>
 * y <- sum(x[])
 * </pre>
 */
class Sum : public ScalarFunc
{
public:
  Sum ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
  bool isLinear(std::vector<bool> const &mask,
                std::vector<bool> const &mask) const;
  bool isScale(unsigned int index, std::vector<bool> const &mask) const;
};

#endif /* FUNC_SUM_H_ */
