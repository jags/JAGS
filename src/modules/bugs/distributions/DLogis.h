#ifndef DLOGIS_H_
#define DLOGIS_H_

#include <distribution/DistReal.h>

/** 
 * Logistic distribution
 * <pre>
 * X ~ dlogis(mu, tau)
 * P(X <= x | mu, tau) = 1/(1 + exp(-(x - mu) * tau))
 * </pre>
 */
class DLogis : public DistReal {
 public:
  DLogis();
  ~DLogis();

  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /**
   * Checks that tau > 0
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DLOGIS_H_ */
