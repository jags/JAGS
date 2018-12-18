#ifndef FUNC_NEG_H_
#define FUNC_NEG_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace base {

/**
 * @short Unary minus operator
 * UnaryMinus implements the prefix operator "-"
 * y <- -x
 */
class Neg : public ScalarFunction
{
public:
    Neg ();
    double evaluate(std::vector<double const *> const &args) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    bool isScale(std::vector<bool> const &mask, 
		 std::vector<bool> const &fix) const;
    std::string deparse(std::vector<std::string> const &par) const;
    bool isDifferentiable(unsigned long i) const;
    double gradient(std::vector<double const *> const &args,
		    unsigned long i) const;

};

}}

#endif /* FUNC_NEG_H_ */
