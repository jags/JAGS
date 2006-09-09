#ifndef DHYPER_H_
#define DHYPER_H_

#include <distribution/DistDiscrete.h>

/**
 * @short Hypergeometric distribution
 * <pre>
 * R ~ dhyper(n1,n2,m1,psi)
 * </pre>
 */
class DHyper : public DistDiscrete {
 public:
  DHyper();
  ~DHyper();
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double x, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  double l(std::vector<SArray const *> const &parameters, bool absolute) const;
  double u(std::vector<SArray const *> const &parameters, bool absolute) const;
  /**
   * Checks that n1, n2, m1 are discrete
   */
  bool checkParameterDiscrete (std::vector<bool> const &mask) const;
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DHYPER_H_ */
