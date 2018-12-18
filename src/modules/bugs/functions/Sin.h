#ifndef FUNC_SIN_H_
#define FUNC_SIN_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Sine function
     * @see Cos
     * <pre>
     * y <- sin(x)
     * </pre>
     */
    class Sin : public ScalarFunction
    {
    public:
	Sin ();
	double evaluate(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_SIN_H_ */
