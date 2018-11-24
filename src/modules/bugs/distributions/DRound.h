#ifndef DROUND_H_
#define DROUND_H_

#include <distribution/ScalarDist.h>

namespace jags {
namespace bugs {

/**
 * @short Rounded observations
 * <pre>
 * x ~ dround(t, d)
 * f(x|t,d) = 1 if x = t rounded to d decimal places
 *          = 0 otherwise
 * </pre>
 */
class DRound : public ScalarDist {
public:
    DRound();
    bool checkParameterValue(std::vector<double const *> const &parameters) 
	const;
    bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    double logDensity(double const x, PDFType type,
		      std::vector<double const *> const &parameters,
		      double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
			double const *lower, double const *upper, 
			RNG *rng) const;
    double l(std::vector<double const *> const &parameters) const;
    double u(std::vector<double const *> const &parameters) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    bool fullRank() const;
    double KL(std::vector<double const *> const &par0,
	      std::vector<double const *> const &par1) const;
};

}}

#endif /* DROUND_H_ */
