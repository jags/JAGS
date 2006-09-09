#ifndef DBERN_H_
#define DBERN_H_

#include <distribution/DistDiscrete.h>

/**
 * @short Bernoulli distribution
 * <pre>
 * R ~ dbern(p)
 * f(r | p) = p^r * (1 - p)^(1 -r) ; r in 0:1
 * </pre>
 */
class DBern : public DistDiscrete {
 public:
  DBern();
  ~DBern();
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double x, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /** Checks that p lies in the open interval (0,1) */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DBERN_H_ */
