#ifndef DMULTI_H_
#define DMULTI_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * <pre>
 * X[] ~ dmulti(p[], N)
 * f(x | p, N) = prod (p^x) ; sum(x) = N
 * </pre>
 * @short Multinomial distribution
 */
class DMulti : public Distribution {
public:
  DMulti();
  ~DMulti();
  double logLikelihood(SArray const &x, 
		       std::vector<SArray const *> const &parameters) const;
  void randomSample(SArray &x, std::vector<SArray const *> const &parameters,
		    RNG *rng) const;
  /**
   * Checks that elements of p lie in range (0,1) and 
   * and sum to 1. Checks that N >= 1
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
  /** Checks that p is a vector and N is a scalar */
  bool checkParameterDim(std::vector<Index> const &dim) const;
  /** Checks that N is discrete-valued */
  bool checkParameterDiscrete(std::vector<bool> const &mask) const;
  Index dim(std::vector<Index> const &dim) const;
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  void typicalValue(SArray &x, std::vector<SArray const *> const &par) const;
};

#endif /* DMULTI_H_ */
