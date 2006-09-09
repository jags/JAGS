#ifndef FUNC_SUBTRACT_H_
#define FUNC_SUBTRACT_H_

#include <function/ScalarFunc.h>

/**
 * @short Difference of two scalars 
 * 
 * The Subtract function implements the infix operator "-" in the
 * BUGS language.  
 * <pre>
 * y <- x - z
 * </pre>
 */
class Subtract : public ScalarFunc 
{ 
public:
  Subtract ();
  double eval (std::vector <SArray const *> const &args) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
  bool isLinear(std::vector<bool> const &mask, 
		std::vector<bool> const &fixmask) const;
};

#endif /* FUNC_SUBTRACT_H_ */
