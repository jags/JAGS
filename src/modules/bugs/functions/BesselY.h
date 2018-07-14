#ifndef BESSEL_Y_H_
#define BESSEL_Y_H_

#include <function/ScalarFunction.h>

namespace jags {
    namespace bugs {
	
	/**
	 * @short Bessel function of the first kind
	 *
	 * <pre>
	 * x <- besselY(x, order)
	 * </pre>
	 **/
	class BesselY : public ScalarFunction
	{
	public:
	    BesselY ();
	    double evaluate(std::vector<double const *> const &args) const;
	    bool checkParameterValue(std::vector <double const *> const &args) const;
	};
	
    }
}

#endif /* BESSEL_Y_H_ */
