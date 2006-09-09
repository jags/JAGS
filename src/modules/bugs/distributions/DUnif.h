#ifndef DUNIF_H_
#define DUNIF_H_

#include <distribution/DistReal.h>

/**
 * <pre>
 * x ~ dunif(a, b)
 * f(x|a,b) = 1/(a - b) ; a <= x <= b
 * </pre>
 * @short Uniform distribution
 */
class DUnif : public DistReal {
 public:
  DUnif();
  ~DUnif();


  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  double l(std::vector<SArray const*> const &parameters, bool absolute) const;
  double u(std::vector<SArray const*> const &parameters, bool absolute) const;
  /** 
   * Checks that a < b
   */
  bool checkParameterValue (std::vector<SArray const *> const &parameters) const;
};

#endif /* DUNIF_H_ */
