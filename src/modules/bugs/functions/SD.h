#ifndef FUNC_SD_H_
#define FUNC_SD_H_

#include <function/ScalarFunc.h>

/**
 * @short Standard deviation function
 * SD calculates the standard deviation of the elements of an array.
 * The array must have at least 2 elements.
 * @see mean
 * <pre>
 * y <- sd(x[])
 * </pre>
 */
class SD : public ScalarFunc
{
public:
  SD ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
};

#endif /* FUNC_SD_H_ */
