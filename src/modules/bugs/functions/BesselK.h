#ifndef BESSEL_K_H_
#define BESSEL_K_H_

#include <function/ScalarFunction.h>

namespace jags {
    namespace bugs {
	
	/**
	 * @short Modified Bessel function of the third kind
	 *
	 * <pre>
	 * x <- besselK(x, order, expon.scaled)
	 * </pre>
	 **/
	class BesselK : public ScalarFunction
	{
	public:
	    BesselK ();
	    double evaluate(std::vector<double const *> const &args) const;
	    bool checkParameterValue(std::vector <double const *> const &args);
	};
	
    }
}

#endif /* BESSEL_K_H_ */
