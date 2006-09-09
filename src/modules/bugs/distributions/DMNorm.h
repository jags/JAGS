#ifndef DMNORM_H_
#define DMNORM_H_

#include <distribution/Distribution.h>

/**
 * @short Multivariate normal distribution
 * <pre>
 * x[] ~ dmnorm(mu[], T[,])
 * f(x | mu, T) = sqrt(det(T)) * exp(-1/2 * (x-mu) %*% T %*% t(x-mu))
 * </pre>
 */
class DMNorm : public Distribution {
public:
  DMNorm();
  ~DMNorm();

  double logLikelihood(SArray const &x,
		       std::vector<SArray const *> const &parameters) const;
  void randomSample(SArray &x,
		    std::vector<SArray const *> const &parameters,
		    RNG *rng) const;
  /**
   * Checks that mu is a vector, T is a square matrix and the sizes of
   * mu and T conform.
   */
  bool checkParameterDim (std::vector<Index> const &dims) const;
  /**
   * Checks that T is symmetric. Note that there is currently NO CHECK
   * that T is positive definite.
   */
  bool checkParameterValue (std::vector<SArray const *> const &parameters) const;
  Index dim(std::vector<Index> const &dims) const;
  /**
   * Convenience random sampler function 
   */
  static void randomsample(double *x, double const *mu, double const *tau,
			   int nrow, RNG *rng);
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  void typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
    const;
};

#endif /* DMNORM_H_ */
