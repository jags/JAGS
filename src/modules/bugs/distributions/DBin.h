#ifndef DBIN_H_
#define DBIN_H_

#include <distribution/DistDiscrete.h>

/**
 * @short Binomial distribution
 * <pre>
 * R ~ dbin(p, n)
 * f(r|p,n) = n!/r!.(n - r)! p^r * (1 - p)^(n - r) ; r in 0:n
 * </pre>
 */
class DBin : public DistDiscrete {
 public:
  DBin();
  ~DBin();
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double x, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  double l(std::vector<SArray const *> const &parameters, bool asbolute) const;
  double u(std::vector<SArray const *> const &parameters, bool absolute) const;
  /**
   * Checks that n is discrete-valued.
   */
  bool checkParameterDiscrete (std::vector<bool> const &mask) const;
  /**
   * Checks that p lies in (0,1) and n > 1
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DBIN_H_ */
