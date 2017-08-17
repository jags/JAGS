#ifndef RE_GAMMA_SLICER2_H_
#define RE_GAMMA_SLICER2_H_

#include <sampler/Slicer.h>

namespace jags {
    namespace glm {

	class REGamma2;

	class REGammaSlicer2 : public Slicer
	{
	    REGamma2 const *_regamma;
	    double const *_shape;
	    double const *_rate;
	    double _sigma, _sigma0;
	  public:
	    REGammaSlicer2(REGamma2 const *regamma,
			   double const *shape, double const *rate,
			   double sigma);
	    double value() const;
	    void setValue(double x);
	    void getLimits(double *lower, double *upper) const;
	    double logDensity() const;
	    void setSigma(double x);
	    void update(RNG *rng);
	};
	
    }
}

#endif /* RE_GAMMA_SLICER2_H_ */
