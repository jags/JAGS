#ifndef DGAMPOIS_H_
#define DGAMPOIS_H_

#include <distribution/RScalarDist.h>

namespace jags {
namespace bugs {

/**
 * <pre>
 * x ~ dgampois(m, s)
 * f(x|m,s) = ((x+s-1)!/(x!*(s-1)!)) * ( s/(s+m) )^s * (1-( s/(s+m)) )^x
 * </pre>
 * @short Gamma-Poisson compound distribution (equivalent to Negative Binomial)
 */
class DGamPois : public RScalarDist {
 public:
  DGamPois();
  double d(double x, PDFType type,
	   std::vector<double const *> const &parameters, 
	   bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that p lies in the interval (0,1) and r > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters) const;
  double KL(std::vector<double const *> const &par0,
	    std::vector<double const *> const &par1) const;
};

}}

#endif /* DGAMPOIS_H_ */
