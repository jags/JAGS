#ifndef DBETA_BIN_H_
#define DBETA_BIN_H_

#include <distribution/RScalarDist.h>

namespace jags {
namespace mix {

/**
 * @short Beta-Binomial distribution
 * <pre>
 * Y ~ dbetabin(a, b, n)
 * f(y|a,b,n) = choose(n, x) * beta(x + a, n - x + b) / beta(a, b);
 * </pre>
 */
class DBetaBin : public RScalarDist {
 public:
  DBetaBin();
  std::string alias() const;
  double d(double x, PDFType type,
	   std::vector<double const *> const &parameters, 
	   bool give_log) const;
  double p(double x, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  double l(std::vector<double const *> const &parameters) const;
  double u(std::vector<double const *> const &parameters) const;
  /**
   * Checks that n is discrete-valued.
   */
  bool checkParameterDiscrete (std::vector<bool> const &mask) const;
  /**
   * Checks that p lies in (0,1),  n >= 1, a > 0, b > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}}

#endif /* DBETA_BIN_H_ */
