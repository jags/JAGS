#ifndef DSAMPLE_H_
#define DSAMPLE_H_

#include <distribution/VectorDist.h>

struct RNG;

namespace jags {
    namespace bugs {

	/**
	 * <pre> X[] ~ dsample(p[], K) </pre> 
	 *
	 * Sample K elements without replacement using a vector of
	 * probability weights p. Sampled values have indicator
	 * 1. Non-sampled values have indicator 0.
	 *
	 * For sampling with replacement see DMulti.
	 *
	 * @short Sample without replacement
	 */
	class DSample : public VectorDist {
	public:
	    DSample();
	    
	    double logDensity(double const *x, unsigned long length,
			      PDFType type, 
			      std::vector<double const *> const &parameters,
			      std::vector<unsigned long> const &lengths,
			      double const *lower, double const *upper) const;
	    void randomSample(double *x, unsigned long length,
			      std::vector<double const *> const &parameters,
			      std::vector<unsigned long> const &lengths,
			      double const *lower, double const *upper,
			      RNG *rng) const;
	    bool checkParameterValue(std::vector<double const *> const &par,
				     std::vector<unsigned long> const &lengths)
		const;
	    bool checkParameterLength(std::vector<unsigned long> const &lengths)
		const;
	    bool checkParameterDiscrete(std::vector<bool> const &mask) const;
	    unsigned long length(std::vector<unsigned long> const &lengths) const;
	    void support(double *lower, double *upper, unsigned long length,
			 std::vector<double const *> const &parameters,
			 std::vector<unsigned long> const &lengths) const;
	    bool isSupportFixed(std::vector<bool> const &fixmask) const;
	    unsigned long df(std::vector<unsigned long> const &lengths) const;
	    bool isDiscreteValued(std::vector<bool> const &mask) const;
	};
	
    }
}

#endif /* DSAMPLE_H_ */
