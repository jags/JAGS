#ifndef SCALAR_FUNCTION_H_
#define SCALAR_FUNCTION_H_

#include <function/Function.h>

namespace jags {

/**
 * @short Scalar-valued Function with scalar arguments
 *
 * ScalarFunction represents scalar-valued functions whose parameters are
 * also scalars. 
 */
class ScalarFunction : public Function
{
public:
    ScalarFunction(std::string const &name, unsigned long npar);
    /**
     * Evaluates the function and returns the result
     *
     * @param args Vector of arguments
     */
    virtual double evaluate(std::vector <double const *> const &args) const = 0;
    /**
     * Checks that the argument values are in the domain of the function.
     */
    virtual bool 
	checkParameterValue(std::vector<double const *> const &args) const;
    /**
     * Tests whether the function preserves power transformations,
     * i.e.  can be expressed as f(x) = a*b^x for some value of a,b,
     * where a,b, and x are all scalar.
     *
     * The default method returns true if the ScalarFunction preserves
     * scale transformations AND only one element of mask is true.
     * Such functions also preserve fixed power transformations.
     *
     * Note that if a > 0 then a power function is a linear function on
     * the log scale: log(f(x)) = log(a) + log(b) * x.
     *
     * @see Function#isPower
     */
    bool isPower(std::vector<bool> const &mask,
		 std::vector<bool> const &isfixed) const;
    /**
     *
     */
    virtual double gradient(std::vector <double const *> const &args,
			    unsigned long i) const;
};

} /* namespace jags */

#endif /* SCALAR_FUNCTION_H_ */
