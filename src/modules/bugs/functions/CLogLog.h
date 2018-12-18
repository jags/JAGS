#ifndef FUNC_CLOGLOG_H_
#define FUNC_CLOGLOG_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short the complementary log log link function
     * @see ICLogLog
     * <pre>
     * y <- cloglog(p)
     * y = log(-log(1-p)) for 0 < p < 1
     * </pre>
     */
    class CLogLog : public ScalarFunction
    {
    public:
	CLogLog ();
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_CLOGLOG_H_ */
