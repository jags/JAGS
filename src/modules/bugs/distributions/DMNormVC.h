#ifndef DMNORM_VC_H_
#define DMNORM_VC_H_

#include <distribution/ArrayDist.h>

namespace jags {

struct RNG;

namespace bugs {

/**
 * @short Multivariate normal distribution with mean-variance parameterization
 * <pre>
 * x[] ~ dmnorm.vcov(mu[], T[,])
 * f(x | mu, T) = exp(-1/2 * (x-mu) %*% inverse(T) %*% t(x-mu)) / sqrt(det(T))
 * </pre>
 */
class DMNormVC : public ArrayDist {
public:
  DMNormVC();

  double logDensity(double const *x, unsigned long length, PDFType type,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims,
		    double const *lower, double const *upper)  const;
  void randomSample(double *x, unsigned long length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims,
		    double const *lower, double const *upper, RNG *rng) const;
  /**
   * Checks that mu is a vector, T is a square matrix and the sizes of
   * mu and T conform.
   */
  bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) 
      const;
  /**
   * Checks that T is symmetric. Note that there is currently NO CHECK
   * that T is positive definite.
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
                           std::vector<std::vector<unsigned long> > const &dims) 
      const;
  std::vector<unsigned long> 
      dim(std::vector<std::vector<unsigned long> > const &dims) const;
  /**
   * Convenience random sampler function 
   *
   * @param x  Array that will hold the result
   *
   * @param mu Array of mean values. A null pointer may be given, and in
   * this case a mean of zero is used.
   *
   * @param tau Either the precision matrix or the variance-covariance
   * matrix.  Only the lower triangle is used.
   * 
   * @param prec Logical flag. If true then tau represents the
   * precision matrix (the parameterization used in the bugs language).
   * If false then tau represents the variance-covariance matrix.
   *
   * @param nrow length of x and mu; number of rows of the square matrix tau
   *
   * @param rng Random number generator
   */
  static void randomsample(double *x, double const *mu, double const *tau,
			   bool prec, int nrow, RNG *rng);
  void support(double *lower, double *upper, unsigned long length,
	       std::vector<double const *> const &parameters,
               std::vector<std::vector<unsigned long> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}}

#endif /* DMNORM_VC_H_ */
