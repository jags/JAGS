#ifndef VECTOR_DIST_H_
#define VECTOR_DIST_H_

#include <distribution/Distribution.h>

#include <vector>
#include <string>

namespace jags {

struct RNG;

/**
 * @short Vector-valued Distribution 
 *
 * The VectorDist class represents distributions for which either the
 * value or the parameters are vectors.
 */
class VectorDist : public Distribution
{
  public:
    /**
     * Constructor.
     * @param name name of the distribution as used in the BUGS language
     * @param npar number of parameters, excluding upper and lower bounds
     */
    VectorDist(std::string const &name, unsigned int npar);
    /**
     * @param x Value at which to evaluate the density.
     *
     * @param length Size of the array x.
     *
     * @param type Indicates whether the full probability density
     * function is required (PDF_FULL) or whether partial calculations
     * are permitted (PDF_PRIOR, PDF_LIKELIHOOD). See PDFType for
     * details.
     * 
     * @param parameters Vector of parameter values of the
     * distribution.
     * 
     * @param lengths Vector of parameter lengths corresponding to the
     * parameter vector.
     *
     * @returns the log probability density.  If the density should be
     * zero because x is inconsistent with the parameters then -Inf is
     * returned. If the parameters are invalid
     * (i.e. checkParameterValue returns false), then the return value
     * is undefined.
     * 
     */
    virtual double 
	logDensity(double const *x, unsigned long length, PDFType type,
		   std::vector<double const *> const &parameters,
		   std::vector<unsigned long> const &lengths,
		   double const *lbound, double const *ubound) const = 0;
    /**
     * Draws a random sample from the distribution. 
     *
     * @param x Array to which the sample values are written
     *
     * @param length Size of the array x.
     *
     * @param parameters  Vector of parameter values at which
     * to evaluate the likelihood. This vector should be of length
     * npar().
     *
     * @param lengths Vector of lengths of the arrays in the argument
     * "parameters".
     * 
     * @param rng pseudo-random number generator to use.
     *
     * @exception length_error 
     */
    virtual void randomSample(double *x, unsigned long length,
			      std::vector<double const *> const &parameters,
			      std::vector<unsigned long> const &lengths, 
			      double const *lbound, double const *ubound,
			      RNG *rng) const = 0;
    /**
     * Returns the support of an unbounded distribution
     */
    virtual void support(double *lower, double *upper, unsigned long length,
			 std::vector<double const *> const &params, 
			 std::vector<unsigned long> const &lengths) const = 0;
    /**
     * Indicates whether the support of the distribution is fixed.
     *
     * @param fixmask Boolean vector of length npar() indicating which
     * parameters have fixed values.
     */
    virtual bool isSupportFixed(std::vector<bool> const &fixmask) const = 0;
    /**
     * Checks that lengths of the parameters are correct.
     */
    virtual bool 
	checkParameterLength (std::vector<unsigned long> const &parameters) 
	const = 0;
    /**
     * Checks that the values of the parameters are consistent with
     * the distribution. For example, some distributions require
     * than certain parameters are positive, or lie in a given
     * range.
     *
     * This function assumes that checkParameterLength returns true.
     */
    virtual bool 
	checkParameterValue(std::vector<double const *> const &parameters,
			    std::vector<unsigned long> const &lengths) const = 0;
    /**
     * Calculates what the length of a sampled value should be, based
     * on the lengths of the parameters.
     *
     * @param par vector of lengths of the parameters.
     */
    virtual unsigned long 
	length (std::vector<unsigned long> const &par) const = 0;
    /**
     * Returns the number of degrees of freedom of the distribution
     * given the parameter lengths. By default this is the same as
     * VectorDist#length. However, some distributions are constrained:
     * and the support occupies a lower dimensional subspace. In this
     * case, the df member function must be overrideen.
     */
    virtual unsigned long df(std::vector<unsigned long> const &lengths) const;
    /**
     * Returns a Monte Carlo estimate of the Kullback-Leibler
     * divergence between distributions with two different parameter
     * values. This is done by drawing random samples from the
     * distribution with the first set of parameters and then
     * calculating the log-likelihood ratio with respect to the second
     * set of parameters.
     *
     * Only one lower and one upper bound is required, which is
     * assumed common to both sets of parameters. This is because the
     * Kullback-Leibler divergence is infinite between two bounded
     * distributions if they do not share the same bounds.
     *
     * @param par1 First set of parameters
     * @param par2 Second set of parameter values
     * @param lengths Vector of parameter lengths, common to both par1 and par2
     * @param lower Pointer to lower bound (NULL if unbounded)
     * @param upper Pointer to upper bound (NULL if unbounded)
     * @param rng Random number generator
     * @param nrep Number of replicates on which to base the estimate
     */
    double KL(std::vector<double const *> const &par1,
	      std::vector<double const *> const &par2,
	      std::vector<unsigned long> const &lengths,
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
  virtual double KL(std::vector<double const *> const &par1,
		    std::vector<double const *> const &par2,
		    std::vector<unsigned long> const &lengths) const;
};

} /* namespace jags */

#endif /* VECTOR_DISTRIBUTION_H_ */

