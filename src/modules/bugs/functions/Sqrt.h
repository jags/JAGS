#ifndef FUNC_SQRT_H_
#define FUNC_SQRT_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Square root function
     * <pre>
     * sigma <- sqrt(tau)
     * tau = sigma * sigma; tau > 0
     * </pre>
     */
    class Sqrt : public ScalarFunction
    {
    public:
	Sqrt ();
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
        bool isPower(std::vector<bool> const &mask,
                     std::vector<bool> const &fix) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_SQRT_H_ */
