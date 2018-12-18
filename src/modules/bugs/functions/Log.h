#ifndef FUNC_LOG_H_
#define FUNC_LOG_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Log function
     * @see Exp
     * <pre>
     * y <- log(x)
     * y = log(x) for x > 0
     * </pre>
     */
    class Log : public ScalarFunction
    {
    public:
	Log ();
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_LOG_H_ */
