#ifndef BESSEL_J_H_
#define BESSEL_J_H_

#include <function/ScalarFunction.h>

namespace jags {
    namespace bugs {
	
	/**
	 * @short Bessel function of the first kind
	 *
	 * <pre>
	 * x <- besselJ(x, order)
	 * </pre>
	 **/
	class BesselJ : public ScalarFunction
	{
	public:
	    BesselJ ();
	    double evaluate(std::vector<double const *> const &args) const;
	    bool checkParameterValue(std::vector <double const *> const &args) const;
	};
	
    }
}

#endif /* BESSEL_J_H_ */
