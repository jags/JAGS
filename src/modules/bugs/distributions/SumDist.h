#ifndef SUM_DIST_H_
#define SUM_DIST_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace bugs {

/**
 * @short Sum of 2 or more random variables
 */
class SumDist : public VectorDist {
public:
    SumDist();

    double logDensity(double const *x, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths) const;
    void randomSample(double *x,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned long> const &lengths,
		      RNG *rng) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    bool fullRank() const;
    bool checkParameterValue(std::vector<double const *> const &params,
			     std::vector<unsigned long> const &lengths) const;
    bool checkParameterLength(std::vector<unsigned long> const &lengths) const;
    bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    void support(double *lower, double *upper, 
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned long> const &lengths) const;
    unsigned long length(std::vector <unsigned long> const &lengths) const;
};

}}

#endif /* SUM_DIST_H_ */
