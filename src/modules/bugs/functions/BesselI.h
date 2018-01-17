#ifndef BESSEL_I_H_
#define BESSEL_I_H_

#include <function/ScalarFunction.h>

namespace jags {
    namespace bugs {
	
	/**
	 * @short Modified Bessel function of the first kind
	 *
	 * <pre>
	 * x <- besselI(x, order, expon.scaled)
	 * </pre>
	 **/
	class BesselI : public ScalarFunction
	{
	public:
	    BesselI ();
	    double evaluate(std::vector<double const *> const &args) const;
	    bool checkParameterValue(std::vector <double const *> const &args);
	};
	
    }
}

#endif /* BESSEL_I_H_ */
