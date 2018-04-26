#ifndef ARRAY_DIST_H_
#define ARRAY_DIST_H_

#include <distribution/Distribution.h>

#include <vector>
#include <string>

namespace jags {

struct RNG;

/**
 * @short Matrix- or array-valued distribution
 *
 * This is the most general sub-class for distributions and is used
 * whenever a distribution takes values in a matrix or
 * higher-dimensional array (e.g. Wishart) or has parameters that are
 * array-valued (e.g. multivariate normal).
 */
class ArrayDist : public Distribution
{
public:
    /**
     * Constructor.
     * @param name name of the distribution as used in the BUGS language
     * @param npar number of parameters, excluding upper and lower bounds
     */
    ArrayDist(std::string const &name, unsigned int npar);
    /**
     * @param x Value at which to evaluate the density.
     *
     * @param type Indicates whether the full probability density
     * function is required (PDF_FULL) or whether partial calculations
     * are permitted (PDF_PRIOR, PDF_LIKELIHOOD). See PDFType for
     * details.
     *
     * @param length Size of the array x.
     *
     * @param parameters Vector of parameter values of the
     * distribution.
     * 
     * @param dims Dimensions of the parameters.
     *
     * @returns The log probability density.  If the density should be
     * zero because x is inconsistent with the parameters then -Inf is
     * returned. If the parameters are invalid
     * (i.e. checkParameterValue returns false), then the return value
     * is undefined.
     */
    virtual double 
	logDensity(double const *x, unsigned long length, PDFType type,
		   std::vector<double const *> const &parameters,
		   std::vector<std::vector<unsigned long> > const &dims,
		   double const *lbound, double const *ubound) const = 0;
    /**
     * Draws a random sample from the distribution. 
     *
     * @param x Array to which the sample values are written
     *
     * @param length Size of the array x.
     * 
     * @param parameters Parameters for the distribution. This vector
     * should be of length npar().  Each element is a pointer to the
     * start of an array containing the parameters. The size of the 
     * array should correspond to the dims parameter. 
     *
     * @param dims Dimensions of the parameters
     *
     * @param lbound pointer to array containing the lower boundary of
     * the distribution. This should be of size length or may be NULL if
     * there is no lower boundary.
     *
     * @param ubound pointer to array containing the upper boundary of
     * the distribution. This should be of size length or may be NULL if
     * there is no upper boundary.
     *
     * @param rng pseudo-random number generator to use.
     *
     * @exception length_error 
     */
    virtual void 
	randomSample(double *x, unsigned long length,
		     std::vector<double const *> const &parameters,
		     std::vector<std::vector<unsigned long> > const  &dims,
		     double const *lbound, double const *ubound,
		     RNG *rng) 	const = 0;
    /**
     * Checks that dimensions of the parameters are correct.
     */
    virtual bool 
	checkParameterDim (std::vector<std::vector<unsigned long> > const &parameters) 
	const = 0;
    /**
     * Checks that the values of the parameters are consistent with
     * the distribution. For example, some distributions require than
     * certain parameters are positive, or lie in a given range.
     *
     * This function assumes that checkParameterDim returns true.
     */
    virtual bool 
	checkParameterValue(std::vector<double const *> const &parameters,
			    std::vector<std::vector<unsigned long> > const &dims) const = 0;
    /**
     * Calculates what the dimension of the distribution should be,
     * based on the dimensions of its parameters. 
     */
    virtual std::vector<unsigned long> 
	dim (std::vector <std::vector<unsigned long> > const &args) const = 0;
    /**
     * Returns the number of degrees of freedom of the distribution
     * given the dimensions of the parameters. By default this is the
     * product of the elements of the dimension vector returned by
     * ArrayDist#dim. However, some distributions are constrained: and
     * the support occupies a lower dimensional subspace. In this
     * case, the df member function must be overrideen.
     */
    virtual unsigned long df(std::vector<std::vector<unsigned long> > const &dims)
	const;
    /**
     * Returns the support of an unbounded distribution
     */
    virtual void 
	support(double *lower, double *upper, unsigned long length,
		std::vector<double const *> const &support,
		std::vector<std::vector<unsigned long> > const &dims) const = 0;
    /**
     * Returns a Monte Carlo estimate of the Kullback-Leibler
     * divergence between distributions with two different parameter
     * values. This is done by drawing random samples from the
     * distribution with the first set of parameters and then
     * calculating the log-likelihood ratio with respect to the second
     * set of parameters.
     *
     * A subclass of ArrayDist can overload this function if the
     * Kullback-Leibler divergence for the distribution it represents
     * can be expressed in closed form.
     *
     * @param par1 First set of parameters
     * @param par2 Second set of parameter values
     * @param dims Vector of parameter dimensions, common to both par1 and par2
     * @param rng Random number generator
     * @param nrep Number of replicates on which to base the estimate
     */
    double KL(std::vector<double const *> const &par1,
	      std::vector<double const *> const &par2,
	      std::vector<std::vector<unsigned long> > const &dims,
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
		      std::vector<std::vector<unsigned long> > const &dims)
	const;
};

} /* namespace jags */

#endif /* ARRAY_DIST_H_ */
