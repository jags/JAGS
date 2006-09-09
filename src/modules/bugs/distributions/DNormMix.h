#ifndef DNORM_MIX_H_
#define DNORM_MIX_H_

#include <distribution/DistReal.h>

/**
 * <pre>
 * x ~ dnorm(mu, tau, p)
 * </pre>
 * @short Mixture of normal distribution
 */
class DNormMix : public DistReal {
 public:
  DNormMix();
  ~DNormMix();

  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  bool checkParameterDim(std::vector<Index> const &dims) const;
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
  void typicalValue(SArray&, std::vector<SArray const *> const &par) const;
};

#endif /* DNORM_MIX_H_ */
