#ifndef DDIRCH_H_
#define DDIRCH_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace bugs {

/**
 * @short Dirichlet distribution
 *
 * Zero shape parameters are allowed.  These represent structural
 * zeros: when x ~ ddirch(alpha) is forward sampled, x[i] = 0 when
 * alpha[i] = 0. 
 *
 * <pre>
 * p[] ~ ddirch(alpha[])
 * f(p | alpha) = C * prod(p^(alpha-1))
 * </pre>
 */
class DDirch : public VectorDist {
public:
    DDirch();
    std::string alias() const;
    double logDensity(double const *x, unsigned int length, PDFType type,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper, RNG *rng) const;
    void typicalValue(double *x, unsigned int length,
		      std::vector<double const *> const &parameter,
		      std::vector<unsigned int> const &lengths,
		      double const *lower, double const *upper) const;
    unsigned int length(std::vector<unsigned int> const &lengths) const;
    /**
     * Checks that alpha is a vector of length at least 2
     */
    bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
    /**
     * Checks that each element of alpha is >= 0.
     *
     * Structural zeros are allowed in the Dirichlet distribution.
     * These are represented by the elements of alpha that are set to
     * zero.  This is permitted only if alpha is fixed and there is at
     * least one non-zero element of alpha.
     */
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<unsigned int> const &lengths)
	const;
    void support(double *lower, double *upper, unsigned int length,
		 std::vector<double const *> const &parameters,
		 std::vector<unsigned int> const &lengths) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned int df(std::vector<unsigned int> const &lengths) const;
    double KL(std::vector<double const *> const &par0,
	      std::vector<double const *> const &par1,
	      std::vector<unsigned int> const &len) const;
};

}}

#endif /* DDIRCH_H_ */
