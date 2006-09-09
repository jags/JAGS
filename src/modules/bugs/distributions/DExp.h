#ifndef DEXP_H_
#define DEXP_H_

#include <distribution/DistReal.h>

/**
 * @short Exponential distribution
 * <pre>
 *  x ~ dexp(lambda)
 *  f(x | lambda) = lambda * exp(-lambda * x) ; x >= 0
 * </pre>
 */
class DExp : public DistReal {
 public:
  DExp();
  ~DExp();

  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /**
   * Checks that lambda > 0
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DEXP_H_ */
