#ifndef FUNC_ARCSIN_H_
#define FUNC_ARCSIN_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverse sine function
     * @see Sin
     * <pre>
     * y <- arcsin(x)
     * </pre>
     */
    class ArcSin : public ScalarFunction
    {
    public:
	ArcSin ();
	std::string alias() const;
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_ARCSIN_H_ */
