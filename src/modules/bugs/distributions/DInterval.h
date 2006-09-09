#ifndef DINTERVAL_H_
#define DINTERVAL_H_

#include <distribution/Distribution.h>

/**
 * @short Interval censored distribution
 * <pre>
 * i ~ dinterval(t, cutpoints[])
 * f(i|t) = 1 if t < cutpoints[i] and t >= cutpoints[i-1]
 *        = 0 otherwise
 * </pre>
 */
class DInterval : public Distribution {
 public:
  DInterval();
  ~DInterval();
  
  double logLikelihood(SArray const &x, 
		       std::vector<SArray const *> const &parameters) const;		
  void randomSample(SArray &x, std::vector<SArray const *> const &parameters,
		    RNG *rng) const;
  /**
   * Checks that t is scalar and that cutpoints is a vector
   */
  bool checkParameterDim(std::vector<Index> const &parameters) const;
  /**
   * Checks that cutpoints are in ascending order
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
  Index dim(std::vector<Index> const &dims) const;
  /**
   * Returns zero
   */
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;

  void typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
    const;
};

#endif /* DINTERVAL_H_ */
