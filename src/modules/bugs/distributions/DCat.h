#ifndef DCAT_H_
#define DCAT_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace bugs {

/**
 * @short Categorical distribution
 * <pre>
 * R ~ dcat(p[])
 * f(r|p[]) = p[r] ; r in 1:length(p)
 * </pre>
 */
class DCat : public VectorDist {
public:
    DCat();

    double logDensity(double const *x, unsigned long length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned long length,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      double const *lbound, double const *ubound,
		      RNG *rng) const;
    /**
     * Checks that all elements of p are positive
     */
    bool checkParameterValue(std::vector<double const*> const &parameters,
			     std::vector<unsigned long> const &lengths) const;
    void support(double *lower, double *upper, unsigned long length,
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned long> const &lengths) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    bool checkParameterLength(std::vector<unsigned long> const &lengths) const;
    unsigned long length(std::vector<unsigned long> const &lengths) const;
    double KL(std::vector<double const *> const &par0,
	      std::vector<double const *> const &par1,
	      std::vector<unsigned long> const &lengths) const;
};

}}

#endif /* DCAT_H_ */
