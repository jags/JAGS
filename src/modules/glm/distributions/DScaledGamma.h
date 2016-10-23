#ifndef DSCALED_GAMMA_H_
#define DSCALED_GAMMA_H_

#include <distribution/RScalarDist.h>

namespace jags {
    namespace glm {
	
	/**
	 * @short Prior distribution for scalar precision
	 *
	 * The scaled F distribution on 1, n degrees of freedom
	 * <pre>
	 * f(x| n1, n2) = s Gamma((n1 + n2)/2) / (Gamma(n1/2) Gamma(n2/2))      
	 *                (n1/n2)^(n1/2) (s x)^(n1/2 - 1)                  
	 *                (1 + (n1/n2) s x)^-(n1 + n2)/2   
	 * </pre>

	 */
	class DScaledGamma : public RScalarDist {
	  public:
	    DScaledGamma();
	    
	    double d(double x, PDFType type,
		     std::vector<double const *> const &parameters, 
		     bool log) const;
	    double p(double x, std::vector<double const *> const &parameters,
		     bool lower, bool log) const;
	    double q(double x, std::vector<double const *> const &parameters,
		     bool lower, bool log) const;
	    double r(std::vector<double const *> const &parameters,
		     RNG *rng) const;
	    /**
	     * Check that s > 0 and n > 0
	     */
	    bool checkParameterValue(std::vector<double const *> const &pars)
		const;
	};

    }
}

#endif /* DSCALED_GAMMA_H_ */
