#ifndef DPAR_H_
#define DPAR_H_

#include <distribution/DistReal.h>

/**
 * <pre>
 * x ~ dpar(alpha, c);
 * f(x|alpha,c) = alpha * c^alpha * x^-(alpha + 1); x > c
 * </pre>
 * @short Pareto distribution
 */
class DPar : public DistReal {
 public:
  DPar();
  ~DPar();

  double d(double x, std::vector<SArray const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<SArray const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<SArray const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<SArray const *> const &parameters, RNG *rng) const;
  double l(std::vector<SArray const*> const &parameters, bool absolute) const;
  double u(std::vector<SArray const*> const &parameters, bool absolute) const;
  /** 
   * Checks that alpha > 0, c > 0
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DPAR_H_ */
