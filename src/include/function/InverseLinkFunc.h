#ifndef INVERSE_LINK_FUNC_H_
#define INVERSE_LINK_FUNC_H_

#include <string>
#include <function/ScalarFunc.h>

/**
 * @short Base class for inverse link functions.
 *
 * An inverse link function is a continuous differentiable scalar function
 * with range (-Inf,Inf).
 *
 * Link functions may be expressed in the BUGS language using an
 * S-style replacement function notation. For example the log link may
 * be expressed as log(y) <- a + b*x or as y <- exp(a + b*x).
 * In this case, the name of the inverse link function is "exp"
 * and its link name is "log".
 */
class InverseLinkFunc : public ScalarFunc
{
    const std::string _link;
public:
    /**
     * Constructor
     * 
     * @param name Function name
     * 
     * @param link Name used by replacement-function notation.
     */
    InverseLinkFunc (std::string const &name, std::string const &link);
    /**
     * Returns the name of the link function
     */
    std::string const &linkName () const;
    /** Link function */
    virtual double link(double mu) const = 0;
    /** Gradient of the inverse link function */
    virtual double grad(double eta) const = 0;
};

#endif
