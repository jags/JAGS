#ifndef DDEXP_H_
#define DDEXP_H_

#include <distribution/DistReal.h>

/**
 * @short Double exponential distribution
 * <pre>
 * x ~ ddexp(mu,tau)
 * f(x|mu,tau) = (tau/2) * exp(-tau * abs(x-mu))
 * </pre>
 */
class DDexp : public DistReal {
public:
  DDexp();
  ~DDexp();
  
  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
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

#endif /* DDEXP_H_ */
