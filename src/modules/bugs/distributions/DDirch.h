#ifndef DDIRCH_H_
#define DDIRCH_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * @short Dirichlet distribution
 *
 * Zero shape parameters are allowed.  These represent structural
 * zeros: when x ~ ddirch(alpha) is forward sampled, x[i] = 0 when
 * alpha[i] = 0. To avoid trapping states in the model, structural
 * zeros are only allowed when the array of shape parameters is
 * fixed.
 *
 * <pre>
 * p[] ~ ddirch(alpha[])
 * f(p | alpha) = C * prod(p^alpha)
 * </pre>
 */
class DDirch : public Distribution {
public:
  DDirch();
  ~DDirch();
  double logLikelihood(SArray const &x,
		       std::vector<SArray const *> const &parameters) const;
  void randomSample(SArray &x, std::vector<SArray const *> const &parameters,
		    RNG *rng) const;
  Index dim(std::vector<Index> const &dims) const;
  /**
   * Checks that alpha is a vector of length at least 2
   */
  bool checkParameterDim(std::vector<Index> const &dims) const;
  /**
   * Checks that each element of alpha is >= 0.
   *
   * Structural zeros are allowed in the Dirichlet distribution.
   * These are represented by the elements of alpha that are set to
   * zero.  This is permitted only if alpha is fixed and there is at
   * least one non-zero element of alpha.
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  void typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
     const;

};

#endif /* DDIRCH_H_ */
