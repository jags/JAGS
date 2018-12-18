#ifndef FUNC_ARCCOSH_H_
#define FUNC_ARCCOSH_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverse hyperbolic cosine function
     * @see cosh
     * <pre>
     * y <- arccosh(x)
     * </pre>
     */
    class ArcCosh : public ScalarFunction
    {
    public:
	ArcCosh ();
	std::string alias() const;
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_ARCCOSH_H_ */
