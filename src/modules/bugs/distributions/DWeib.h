#ifndef DWEIB_H_
#define DWEIB_H_

#include <distribution/DistReal.h>

/**
 * <pre>
 * x ~ dweib(a, b)
 * f(x|b,a) = a * b * x^(b - 1) * exp (- a * x^b)
 * </pre>
 * @short Weibull distribution
 */
class DWeib : public DistReal {
public:
  DWeib();
  ~DWeib();

  double d(double x, std::vector<SArray const *> const &parameters, 
	   bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  /** 
   * Checks that a > 0, b > 0
   */
  bool checkParameterValue (std::vector<SArray const *> const &parameters) const;
};

#endif /* DWEIB_H_ */
