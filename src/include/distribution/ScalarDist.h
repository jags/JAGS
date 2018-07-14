#ifndef SCALAR_DIST_H_
#define SCALAR_DIST_H_

#include <distribution/Distribution.h>

namespace jags {

struct RNG;

/**
 * Enumerates three possible ranges of support for a scalar random
 * variable.
 *
 * DIST_UNBOUNDED for support on the whole real line
 *
 * DIST_POSITIVE for support on values > 0
 *
 * DIST_PROPORTION for support on values in [0,1]
 *
 * DIST_SPECIAL for other distributions, e.g. distributions where the
 * support depends on the parameter values.
 */
enum Support {DIST_UNBOUNDED, DIST_POSITIVE, DIST_PROPORTION, DIST_SPECIAL};

/**
 * @short Scalar distributions
 *
 * The ScalarDist class represents distributions that are scalar-valued
 * and for which the parameters are also scalars.
 */
class ScalarDist : public Distribution
{
  const Support _support;
 public:
  /**
   * Constructor
   *
   * @param name BUGS language name of distribution
   *
   * @param npar Number of parameters, excluding upper and lower bound
   *
   * @param support Support of distribution
   */
  ScalarDist(std::string const &name, unsigned int npar, Support support);
  /**
   * This implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  /**
   * Lower limit of distribution, given parameters.  If the
   * distribution has no lower limit, this should return JAGS_NEGINF.
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double l(std::vector<double const *> const &parameters) const;
  /**
   * Upper limit of distribution, given parameters. If the
   * distribution has no upper limit, this should return JAGS_POSINF.
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double u(std::vector<double const *> const &parameters) const;
  /**
   * Calculates log probability density of the distribution
   */
  virtual double logDensity(double x, PDFType type,
			    std::vector<double const *> const &parameters,
			    double const *lbound, double const *ubound)
      const = 0;
  /**
   * Draws a random sample 
   */
  virtual double randomSample(std::vector<double const *> const &parameters, 
			      double const *lbound, double const *ubound,
			      RNG *rng) const = 0;
  /**
   * Checks whether the parameter values are valid
   */
  virtual bool checkParameterValue(std::vector<double const *> const &params) 
      const = 0;
  /**
   * Returns the number of degrees of freedom of the distribution.
   * By default this is 1. For scalar distributions that are
   * deterministic functions of the parameters, this must be overridden.
   */
  virtual unsigned long df() const;
  /**
   * Returns a Monte Carlo estimate of the Kullback-Leibler divergence
   * between distributions with two different parameter values. This
   * is done by drawing random samples from the distribution with the
   * first set of parameters and then calculating the log-likelihood
   * ratio with respect to the second set of parameters.
   *
   * Only one lower and one upper bound is required, which is assumed
   * common to both sets of parameters. This is because the
   * Kullback-Leibler divergence is infinite between two bounded
   * distributions if they do not share the same bounds.
   *
   * @param par1 First set of parameters
   * @param par2 Second set of parameter values
   * @param lower Pointer to lower bound (NULL if there is no lower bound)
   * @param upper Pointer to upper bound (NULL if there is no upper bound)
   * @param rng Random number generator
   * @param nrep Number of replicates on which to base the estimate
   */
  double KL(std::vector<double const *> const &par1,
	    std::vector<double const *> const &par2,
	    double const *lower, double const *upper,
	    RNG *rng, unsigned int nrep) const;
  /**
   * Returns the Kullback-Leibler divergence between distributions
   * with two different parameter values.
   *
   * This is a virtual function that must be overloaded for any
   * distribution that allows exact Kullback-Leibler divergence
   * calculations. The default method returns JAGS_NA, indicating that
   * the method is not implemented.
   */
  virtual double KL(std::vector<double const *> const &,
		    std::vector<double const *> const &) const;
  
};

} /* namespace jags */

#endif /* SCALAR_DIST_H_ */

