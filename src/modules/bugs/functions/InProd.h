#ifndef FUNC_INPROD_H_
#define FUNC_INPROD_H_

#include <function/ScalarFunc.h>

/**
 * @short Inner product of two vectors
 * The InProd function returns the sum of the product of corresponding
 * elements from two arrays with the same dimension
 * <pre>
 * y <- inprod(a[],b[])
 * </pre>
 */
class InProd : public ScalarFunc
{
public:
  InProd ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
  bool isLinear(std::vector<bool> const &mask,
                std::vector<bool> const &fix) const;
  bool isScale(unsigned int index, std::vector<bool> const &fix) const;
};

#endif /* FUNC_INPROD_H_ */
