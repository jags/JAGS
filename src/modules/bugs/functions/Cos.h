#ifndef FUNC_COS_H_
#define FUNC_COS_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Cosine function
     * @see Sin
     * <pre>
     * y <- cos(x)
     * </pre>
     */
    class Cos : public ScalarFunction
    {
    public:
	Cos ();
	double evaluate(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_COS_H_ */
