#ifndef DINTERVAL_H_
#define DINTERVAL_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace bugs {

/**
 * @short Interval censored distribution
 * <pre>
 * i ~ dinterval(t, cutpoints[])
 * f(i|t) = 1 if t < cutpoints[i] and t >= cutpoints[i-1]
 *        = 0 otherwise
 * </pre>
 */
class DInterval : public VectorDist {
public:
    DInterval();
  
    double logDensity(double const *x, unsigned long length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned long length,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lower, double const *upper, RNG *rng) const;
    /**
     * Checks that cutpoints are in ascending order
     */
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<unsigned long> const &lengths) const;
    bool checkParameterLength(std::vector<unsigned long> const &lengths) const;
    void support(double *lower, double *upper, unsigned long length,
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned long> const &lengths) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned long df(std::vector<unsigned long> const &lengths) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    unsigned long length(std::vector<unsigned long> const &params) const;
    double KL(std::vector<double const *> const &par0,
	      std::vector<double const *> const &par1,
	      std::vector<unsigned long> const &lengths) const;
};

}}

#endif /* DINTERVAL_H_ */
