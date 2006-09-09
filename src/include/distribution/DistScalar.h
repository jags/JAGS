#ifndef DIST_SCALAR_H_
#define DIST_SCALAR_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * Enumerates three possible ranges of support for a real valued
 * random variable
 *
 * DIST_UNBOUNDED for support on the whole real line
 *
 * DIST_POSITIVE for support on values > 0
 *
 * DIST_PROPORTION for support on values in [0,1]
 *
 * DIST_SPECIAL for other distributions, i.e. distributions where the
 * support depends on the parameter values.
 */
enum Support {DIST_UNBOUNDED, DIST_POSITIVE, DIST_PROPORTION, DIST_SPECIAL};

/**
 * Abstract class for scalar valued distributions.
 *
 * A subclass of DistScalar has to implement the d,p,q, and r virtual
 * member functions. These are based on the d-p-q-r functions provided
 * by libRmath and, for most of the standard JAGS distributions are
 * simply wrappers around them.
 *
 * The JAGS versions of most (but not all) scalar distributions extend
 * the distribution families in libRmath by allowing the distribution
 * to be bounded.
 *
 * @short Real valued distributions
 */
class DistScalar : public Distribution
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
   *
   * @param boundable Logical flag indicating whether the distribution
   * can be bounded using the T(,) construct.
   *
   * @param discrete Logical flag indicating whether the distribution
   * is discrete-valued.
   */
  DistScalar(std::string const &name, unsigned int npar,
	     Support support, bool canbound, bool discrete);
  /**
   * All scalar distributions have the same dimension
   */
  Index dim(std::vector<Index> const &parameters) const;
  double lowerSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  double upperSupport(unsigned long i,
		      std::vector<SArray const *> const &parameters) const;
  unsigned long df(std::vector<SArray const *> const &parameters) const;
  /**
   * Lower limit of distribution, given parameters, but ignoring
   * bounds.  If the distribution has no lower limit, this should
   * return -DBL_MAX.
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double l(std::vector<SArray const *> const &parameters,
                   bool absolute) const;
  /**
   * Upper limit of distribution, given parameters, but ignoring
   * bounds.  If the distribution has no upper limit, this should
   * return DBL_MAX
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double u(std::vector<SArray const *> const &parameters,
                   bool absolute) const;
  /** 
   * Density function, ignoring bounds
   * @param x value at which to evaluate the density
   * @param parameters Array of parameters
   * @param give_log Indicates whether to return log density. 
   */
  virtual double d(double x, std::vector<SArray const *> const &parameters, 
		   bool give_log) const = 0;
  /** 
   * Distribution function, ignoring bounds
   * @param x quantile at which to evaluate the distribution function
   * @param parameters Array of parameters
   * @param lower If true, return value is P[X <= x]. Otherwise
   * P[X > x]
   * @param give_log Indicates whether to return log probabability
   */
  virtual double p(double x, std::vector<SArray const *> const &parameters, bool lower,
		   bool give_log) const = 0;
  /** 
   * Quantile function, ignoring bounds
   * @param p probability for which to evaluate quantile
   * @param parameters Array of parameters
   * @param log_p Indicates whether p is given as log(p). 
   */
  virtual double q(double p, std::vector<SArray const *> const &parameters, 
		   bool lower, bool log_p) const = 0;
  /**
   * Random number generation, ignoring bounds
   * @param parameters Array of parameters
   */
  virtual double 
    r(std::vector<SArray const *> const &parameters, RNG *rng) const = 0;
  /**
   * Returns the median. Note that this function can be overloaded
   * by a subclass if necessary.
   */
  void typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
    const;
  /**
   * Most scalar-valued distributions have only scalar parameters. For
   * convenience, a default checkParameterDim member function is provided
   * that checks returns true if the value and the parameters are all
   * scalar.  Distributions that do not take scalar parameters need to
   * overload this
   */
  bool checkParameterDim (std::vector<Index> const &dims) const;
};

#endif /* DIST_SCALAR_H_ */
