#ifndef FUNC_ROUND_H_
#define FUNC_ROUND_H_

#include <function/ScalarFunc.h>

namespace bugs {

    /**
     * @short Round to integer, away from zero
     * <pre>
     * y <- round(x)
     * </pre>
     */
    class Round : public ScalarFunc
    {
    public:
	Round ();
	double evaluateScalar(std::vector<double const *> const &args) const;
	bool isDiscreteValued(std::vector<bool> const &mask) const;
    };
}

#endif /* FUNC_ROUND_H_ */
