#ifndef FUNC_LOGGAM_H_
#define FUNC_LOGGAM_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Log gamma function
     * <pre>
     * y <- loggam(x)
     * y = log(gamma(x)) for x > 0
     * </pre>
     */
    class LogGam:public ScalarFunction
    {
    public:
	LogGam ();
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_LOGGAM_H_ */
