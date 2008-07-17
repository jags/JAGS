#ifndef FUNC_EXP_H_
#define FUNC_EXP_H_

#include <function/InverseLinkFunc.h>

namespace bugs {

    /**
     * @short Exponential function
     * @see Log
     * <pre>
     * y <- exp(x)
     * </pre>
     */
    class Exp: public InverseLinkFunc
    {
    public:
	Exp ();
	double evaluateScalar(std::vector<double const *> const &args) const;
	double link(double mu) const;
	double grad(double eta) const;
    };

}

#endif /* FUNC_EXP_H_ */
