#ifndef FUNC_TANH_H_
#define FUNC_TANH_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Hyperbolic tangent function
     * @see tanh
     * <pre>
     * y <- tanh(x)
     * </pre>
     */
    class Tanh : public ScalarFunction
    {
    public:
	Tanh ();
	double evaluate(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	bool gradient(double &grad, std::vector<double const *> const &args,
		      unsigned long i) const;
    };

}}

#endif /* FUNC_TANH_H_ */
