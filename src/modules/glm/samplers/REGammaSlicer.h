#ifndef RE_GAMMA_SLICER_H_
#define RE_GAMMA_SLICER_H_

#include <sampler/Slicer.h>

#include <vector>

extern "C" {
#include <cholmod.h>
}

using std::vector;

namespace jags {
    namespace glm {

	class Outcome;

	class REGammaSlicer : public Slicer
	{
	    std::vector<Outcome const *> _outcomes;
	    cholmod_sparse const *_x;
	    cholmod_dense const *_z;
	    double const *_shape;
	    double const *_rate;
	    double _sigma, _sigma0;
	    REGammaSlicer *_slicer;
	  public:
	    REGammaSlicer(std::vector<Outcome*> const &outcomes,
			  cholmod_sparse const *x, cholmod_dense const *z,
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

#endif /* RE_GAMMA_SLICER_H_ */
