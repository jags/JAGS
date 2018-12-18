#ifndef FUNC_SINH_H_
#define FUNC_SINH_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Hyperbolic sine function
     * @see sinh
     * <pre>
     * y <- sinh(x)
     * </pre>
     */
    class Sinh : public ScalarFunction
    {
    public:
	Sinh ();
	double evaluate(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_SINH_H_ */
