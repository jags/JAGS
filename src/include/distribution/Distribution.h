#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <vector>
#include <string>

namespace jags {

struct RNG;

/**
 * Specifies different ways of calculating the log probability density
 * function. Typically we are only interested in calculating ratios of
 * probability density functions. When some terms in the log density
 * are constant (and therefore cancel out when taking ratios) they may
 * optionally be omitted from the calculations for improved
 * efficiency.
 *
 * PDF_FULL is used when neither the sampled values nor the
 * pararameters are constant. This implies full evaluation of the log
 * density.
 *
 * PDF_PRIOR is used when the parameters are constant. Terms that
 * depend only on the parameters may be omitted from the log density.
 * 
 * PDF_LIKELIHOOD is used when the sampled value is constant. Terms
 * that depend only on the sampled value may be omitted from the log
 * density.
 *
 * Note that any function taking a PDFType argument is free to ignore
 * it and can always calculate the full density function.
 */
enum PDFType {PDF_FULL, PDF_PRIOR, PDF_LIKELIHOOD};

/**
 * @short Distribution of a random variable
 *
 * Distribution objects contain only constant data members and all
 * member functions are constant. Hence only one object needs to be
 * instantiated for each subclass.
 *
 * The DistTab class provides a convenient way of storing Distribution
 * objects and referencing them by name.
 *
 * @see DistTab
 */
class Distribution
{
    const std::string _name;
    const unsigned long _npar;
public:
    /**
     * Constructor.
     * @param name name of the distribution as used in the BUGS language
     * @param npar number of parameters, excluding upper and lower bounds
     */
    Distribution(std::string const &name, unsigned long npar);
    virtual ~Distribution();
    /**
     * @returns the BUGS language name of the distribution
     */
    std::string const &name() const;
    /**
     * Returns an alternate name for the distribution. The default
     * implementation returns an empty string. However, this function
     * may be overridden to return a non-empty string giving an
     * alternate name for the distribution.
     */
    virtual std::string alias() const;
    /**
     * Indicates whether the support of the distribution is fixed.
     *
     * @param fixmask Boolean vector indicating which parameters have
     * fixed values.
     */
    virtual bool isSupportFixed(std::vector<bool> const &fixmask) const = 0;
    /**
     * Returns the number of parameters of the distribution
     */
    unsigned long npar() const;
    /**
     * Some distributions require some of the parameters to be discrete
     * valued. As most distributions do not require discrete valued paremeters,
     * a default implementation is provided which always returns true.
     *
     * @param mask Boolean vector indicating which parameters are
     * discrete valued.
     */
    virtual bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    /**
     * Returns true if the distribution has support on the integers.The
     * default implementation returns false, so this must be overridden
     * for discrete-valued distributions.
     *
     * @param mask Vector indicating whether parameters are discrete
     * or not. Most implementations will ignore this argument, as a
     * distribution normally has support either on the real line or on
     * the integers. However, this argument is required in order to
     * support observable functions, for which the support may depend
     * on the arguments.
     *
     * @see Function#isDiscreteValued
     */
    virtual bool isDiscreteValued(std::vector<bool> const &mask) const;
    /**
     * Tests for a location parameter.  A parameter of a distribution
     * is considered to be a location parameter if, when it's value is
     * incremented by X, the whole distribution is shifted by X,
     * indpendently of the other parameter values.
     * 
     * This is a virtual function, for which the default implementation
     * always returns false. Distributions with location parameters must
     * overload this function.
     *
     * @param index Index number (starting from 0) of the parameter to be
     * tested.
     */
    virtual bool isLocationParameter(unsigned int index) const;
    /**
     * Tests for a scale parameter.  A parameter of a distribution is
     * considered to be a scale parameter if, when it's value is
     * multiplied by X, the whole distribution multiplied by X,
     * indpendently of the other parameter values.
     * 
     * Note that this definition excludes "location-scale" models:
     * i.e. if the density of y takes the form (1/b)*f((y-a)/b) then b
     * is not considered a scale parameter.
     *
     * This is a virtual function, for which the default
     * implementation always returns false. Distributions with scale
     * parameters must overload this function.
     *
     * @param index Index number (starting from 0) of the parameter to
     * be tested.
     */
    virtual bool isScaleParameter(unsigned int index) const;
    /**
     * Indicates whether the distribution can be bounded. The default
     * implementation returns false.
     */
    virtual bool canBound() const;
};

/**
 * Checks that the distribution accepts npar parameters
 */
inline bool checkNPar(Distribution const *dist, unsigned long npar)
{
    return (dist->npar() == 0 && npar > 0) ||  dist->npar() == npar;
}

} /* namespace jags */

#endif /* DISTRIBUTION_H_ */
