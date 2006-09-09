#ifndef INVERSE_H_
#define INVERSE_H_

#include <function/Function.h>

/**
 * @short Inverts a square matrix
 * <pre>
 * y[,] <- inverse(x[,])
 * </pre>
 */
class Inverse:public Function
{
public:
  Inverse ();

  void evaluate (SArray & value,
		 std::vector <SArray const *> const &args) const;
  Index dim (std::vector <Index> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
};

#endif /* INVERSE_H_ */
