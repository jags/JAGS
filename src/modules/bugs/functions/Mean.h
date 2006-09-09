#ifndef FUNC_MEAN_H_
#define FUNC_MEAN_H_

#include <function/ScalarFunc.h>

/**
 * @short Mean function
 * Mean calculates the mean of the elements of an array
 * @see SD
 * <pre>
 * y <- mean(x[])
 * </pre>
 */
class Mean : public ScalarFunc
{
public:
  Mean ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  bool isLinear(std::vector<bool> const &mask,
                std::vector<bool> const &fix) const;
  bool isScale(unsigned int index, std::vector<bool> const &fix) const;
};

#endif /* FUNC_MEAN_H_ */
