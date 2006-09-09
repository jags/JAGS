#ifndef DCHISQ_H_
#define DCHISQ_H_

#include <distribution/DistReal.h>

/**
 * @short Chi square distribution
 * <pre>
 * x ~ dchisq(k)
 * f(x|k) = 2^(-k/2) * x^(k/2 - 1) * exp(-x/2) / gamma(x/2); k > 0
 * </pre>
 */
class DChisqr : public DistReal {
 public:
  DChisqr();
  ~DChisqr();
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /**
   * Checks that k > 0
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DCHISQ_H_ */
