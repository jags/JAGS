#ifndef DPOIS_H_
#define DPOIS_H_

#include <distribution/DistDiscrete.h>

/**
 * <pre>
 * r ~ dpois(lambda)
 * f(r|lambda) = exp(-lambda) * lambda^r / r! ; r = 0, 1 ,...
 * </pre>
 * @short Poisson distribution
 */
class DPois : public DistDiscrete {
 public:
  DPois();
  ~DPois();

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

#endif /* DPOIS_H_ */
