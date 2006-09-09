#ifndef FUNC_NEG_H_
#define FUNC_NEG_H_

#include <function/ScalarFunc.h>

/**
 * @short Unary minus operator
 * UnaryMinus implements the prefix operator "-"
 * y <- -x
 */
class Neg : public ScalarFunc 
{
public:
  Neg ();
  double eval (std::vector <SArray const *> const &args) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
  bool isLinear(std::vector<bool> const &mask,
		std::vector<bool> const &fixmask) const;
  bool isScale(unsigned int index, std::vector<bool> const &fixmask) const;
};

#endif /* FUNC_NEG_H_ */
