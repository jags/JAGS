#ifndef FUNC_ARCCOS_H_
#define FUNC_ARCCOS_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace bugs {

    /**
     * @short Inverse cosine function
     * @see Cos
     * <pre>
     * y <- arccos(x)
     * </pre>
     */
    class ArcCos : public ScalarFunction
    {
    public:
	ArcCos ();
	std::string alias() const;
	double evaluate(std::vector<double const *> const &args) const;
	bool checkParameterValue(std::vector<double const *> const &args) const;
	bool isDifferentiable(unsigned long i);
	double gradient(std::vector<double const *> const &args,
			unsigned long i) const;
    };

}}

#endif /* FUNC_ARCCOS_H_ */
