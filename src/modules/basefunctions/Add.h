#ifndef FUNC_ADD_H_
#define FUNC_ADD_H_

#include <function/ScalarFunc.h>

/**
 * The Add function implements the infix operator "+" in the BUGS language.
 * The expression a + b + c is broken down into ((a + b) + c), so two
 * separate calls to Add are made.
 * @short Sum of two scalars
 * <pre>
 * y <- a + b
 * </pre>
 */
class Add:public ScalarFunc
{
public:
  Add ();
  double eval (std::vector < SArray const *>const &args) const;
  bool isDiscreteValued(std::vector<bool> const &flags) const;
  bool isLinear(std::vector<bool> const &mask,
		std::vector<bool> const &fixmask) const;
};

#endif /* FUNC_ADD_H_ */
