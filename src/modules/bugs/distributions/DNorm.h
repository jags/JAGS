#ifndef DNORM_H_
#define DNORM_H_

#include <distribution/DistReal.h>

/**
 * <pre>
 * x ~ dnorm(mu, tau)
 * f(x | mu, tau) = sqrt(tau) * exp(-1/2 * tau * (x - mu)^2)
 * </pre>
 * @short Normal distribution
 */
class DNorm : public DistReal {
 public:
  DNorm();
  ~DNorm();

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

#endif /* DNORM_H_ */
