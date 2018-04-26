#ifndef DNORM_MIX_H_
#define DNORM_MIX_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace mix {

/**
 * @short Mixture of normal distribution
 * <pre>
 * x ~ dnorm(mu, tau, pi)
 * </pre>
 * The parameters mu, tau, pi are conforming vectors or arrays. Each
 * element corresponds to a component of the mixture. If all parameters
 * are vectors then component i is normal with mean mu[i], tau[i] and 
 * probability of selection pi[i]/sum(pi).
 *
 * Note that the parameter pi is redundant: the likelihood depends
 * only on pi/sum(pi).
 */
    class DNormMix : public VectorDist {
    public:
	DNormMix();
	
	double logDensity(double const *x, unsigned long length, PDFType type,
			  std::vector<double const *> const &parameters,
			  std::vector<unsigned long> const &lengths,
			  double const *lower, double const *upper) const;
	void randomSample(double *x, unsigned long length,
			  std::vector<double const *> const &parameters,
			  std::vector<unsigned long> const &lengths,
			  double const *lower, double const *upper, RNG *rng) 
	    const;
	void support(double *lower, double *upper, unsigned long length,
		     std::vector<double const *> const &parameters,
		     std::vector<unsigned long> const &lengths) const;
	bool isSupportFixed(std::vector<bool> const &fixmask) const;
	bool checkParameterLength(std::vector<unsigned long> const &lengths)
	    const;
	bool checkParameterValue(std::vector<double const *> const &parameters,
				 std::vector<unsigned long> const &lengths) 
	    const;
	void typicalValue(double *x, std::vector<double const *> const &par,
			  std::vector<unsigned long> const &lengths)  const;
	unsigned long length(std::vector<unsigned long> const &parlengths) const;
    };

}}


#endif /* DNORM_MIX_H_ */
