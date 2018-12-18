#ifndef FUNC_ARCTAN_H_
#define FUNC_ARCTAN_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverse tangent function
     * @see Tan
     * <pre>
     * y <- arctan(x)
     * </pre>
     */
    class ArcTan : public ScalarFunction
    {
    public:
	ArcTan ();
	std::string alias() const;
	double evaluate(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i) const;
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_ARCTAN_H_ */
