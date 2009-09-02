#ifndef FUNC_TRUNC_H_
#define FUNC_TRUNC_H_

#include <function/ScalarFunc.h>

namespace bugs {

/**
 * @short Round to integer, towards zero
 * <pre>
 * y <- trunc(x)
 * </pre>
 */
    class Trunc : public ScalarFunc
    {
    public:
	Trunc ();
	double evaluateScalar(std::vector<double const *> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

}

#endif /* FUNC_TRUNC_H_ */
