#ifndef DGAMMA_H_
#define DGAMMA_H_

#include <distribution/DistReal.h>

/**
 * @short gamma distribution
 * <pre>
 * X ~ dgamma(r, mu)
 * f(x|r,mu) = mu^r * x^(r - 1) * exp(-mu * x)
 * </pre>
 */
class DGamma : public DistReal {
 public:
  DGamma();
  ~DGamma();

  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /**
   * Checks that r > 0, mu > 0
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
  /**
   * The mean
   */
  void typicalValue(SArray &x, std::vector<SArray const *> const &par) const;
};

#endif /* DGAMMA_H_ */
