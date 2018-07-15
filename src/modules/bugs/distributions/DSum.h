#ifndef DSUM_H_
#define DSUM_H_

#include <distribution/ArrayDist.h>

namespace jags {
namespace bugs {

/**
 * @short Sum of 2 or more random variables
 */
class DSum : public ArrayDist {
public:
    DSum();

    double logDensity(double const *x, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned long> > const &dims)
	const;
    void randomSample(double *x,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned long> > const &dims,
		      RNG *rng) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    unsigned long df(std::vector<std::vector<unsigned long> > const &dims) const;
    bool checkParameterValue(std::vector<double const *> const &params,
			     std::vector<std::vector<unsigned long> > const &dims) const;
    bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims)
	const;
    bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    void support(double *lower, double *upper,
		 std::vector<double const *> const &parameters,
		 std::vector<std::vector<unsigned long> > const &dims) const;
    std::vector<unsigned long> 
	dim(std::vector <std::vector<unsigned long> > const &dims) const;
};

}}

#endif /* DSUM_H_ */
