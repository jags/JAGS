#ifndef DCAT_H_
#define DCAT_H_

#include <distribution/DistDiscrete.h>

/**
 * @short Categorical distribution
 * <pre>
 * R ~ dcat(p[])
 * f(r|p[]) = p[r] ; r in 1:dim(p)
 * </pre>
 */
class DCat : public DistDiscrete {
public:
  DCat();
  ~DCat();
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
   * Checks that  p is a vector of length at least 2
   */
  bool checkParameterDim(std::vector<Index> const &dims) const;
  /**
   * Checks that all elements of p lie in the interval (0,1)
   * and sum to 1
   */
  bool checkParameterValue(std::vector<SArray const *> const &parameters) const;
};

#endif /* DCAT_H_ */
