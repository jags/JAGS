#ifndef DMSTATE_H_
#define DMSTATE_H_

#include <distribution/ArrayDist.h>

namespace jags {
namespace msm {

/**
 * Distribution of x in multi-state homogeneous Markov model at time t
 * when state at time 0 was xold and transition intensity matrix is
 * Lambda
 * <pre> 
 * x ~ dmstate(xold, dt, Lambda[,]) 
 * </pre>
 * @short Interval censored transitions in a multistate Markov model
 */
class DMState : public ArrayDist {
public:
    DMState();
  
    double logDensity(double const *x, unsigned long length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned long> > const &dims,
		      double const *lower, double const *upper)	const;
    void randomSample(double *x, unsigned long length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned long> > const &dims,
		      double const *lower, double const *upper, RNG *rng) const;
    void support(double *lower, double *upper, unsigned long length,
		 std::vector<double const *> const &parameters,
		 std::vector<std::vector<unsigned long> > const &dims) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims)
	const;
    bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    bool checkParameterValue(std::vector<double const *> const &par,
			     std::vector<std::vector<unsigned long> > const &dims) const;
    std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> >
				  const &dims) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
};

}}

#endif /* DMSTATE_H_ */

