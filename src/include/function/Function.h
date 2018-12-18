#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <string>
#include <vector>

namespace jags {

class Range;

/**
 * Base class for functions.
 *
 * Function objects contain only constant data members and all member
 * functions are constant. Hence only one object needs to be
 * instantiated for each subclass.
 *
 * The FuncTab class provides a convenient way of storing Function
 * objects and referencing them by name.
 *
 * @see FuncTab
 * @short BUGS-language function or operator
 */
class Function
{
    const std::string _name;
    const unsigned long _npar;
public:
    /**
     * Constructor.
     * @param name name of the function as used in the BUGS language
     * @param npar number of parameters. If npar == 0 then the function
     * takes a variable number of parameters.
     */
    Function(std::string const &name, unsigned long npar);
    virtual ~Function ();
    /**
     * Returns the  BUGS language name of the function
     */
    std::string const &name () const;
   /**
     * Returns an alternate name for the function. The default
     * implementation returns an empty string. However, this function
     * may be overridden to return a non-empty string giving an
     * alternate name for the function.
     */
    virtual std::string alias() const;
    /**
     * Returns the number of parameters. A variadic function returns 0. 
     */
    unsigned long npar() const;
    /**
     * Returns true if the function returns integer values. The
     * default implementation returns false. A function that returns
     * integer values will need to overload this.
     *
     * @param mask Vector indicating whether parameters are discrete
     * or not.
     */
    virtual bool isDiscreteValued(std::vector<bool> const &mask) const;
    /**
     * Checks whether the required arguments are discrete-valued.
     * Most functions do not require discrete-valued arguments, so
     * a default implementation is provided which always returns true.
     *
     * @param mask Boolean vector indicating which parameters are
     * discrete valued.
     */
    virtual bool checkParameterDiscrete(std::vector<bool> const &mask) const;
    /**
     * Tests whether the function preserves scale transformations, i.e.
     * it can be represented element-wise in the form f(x)[i] = B[i] * x[i].
     *
     * It is permitted for the isScale function to return false
     * negatives, i.e. to return false when the function does preserve
     * scale transformations.  Since most functions are not scale
     * transformations, the default method returns false.
     * 
     * @param mask boolean vector of length equal to the number of
     * parameters.  The mask indicates a subset of parameters (those
     * with value true) that are themselves scale transformations of
     * x.  The others are not functions of x. At least one element of
     * mask must be true.
     *
     * @param isfixed boolean vector. This may be empty, in which case
     * it is ignored.  A non-empty vector must have length equal to
     * mask, and denotes the parameters whose values are fixed. In
     * this case the test is more restrictive, and returns true only
     * if the coefficient B[i] is fixed for all i, assuming that the
     * parameters for which mask is true are themselves fixed scale
     * transformations.
     *
     * @see DeterministicNode#isClosed
     */
    virtual bool isScale(std::vector<bool> const &mask,
			 std::vector<bool> const &isfixed) const;
    /**
     * Tests whether the function preserves linear transformations, i.e. 
     * a function of the form f(x) = A + B %*% x.
     *
     * All scale functions are linear functions. Therefore, the
     * default method calls isScale. It only needs to be overridden
     * for linear functions that are not scale functions.
     *
     * @see Function#isScale
     */
    virtual bool isLinear(std::vector<bool> const &mask,
			  std::vector<bool> const &isfixed) const;
    /**
     * Tests whether the function preserves power transformations,
     * i.e.  can be expressed as f(x) = a*b^x for some value of a,b,
     * where a,b, and x are all scalar.
     *
     * The default method returns false.
     */
    virtual bool isPower(std::vector<bool> const &mask,
			 std::vector<bool> const &isfixed) const;
    /**
     * Tests whether the function is additive, i.e.  a function of the
     * form f(x) = A + sum(x) where sum(x) is the sum of all elements of x.
     *
     * The default method returns false
     */
    virtual bool isAdditive(std::vector<bool> const &mask,
			    std::vector<bool> const &isfixed) const;
    /**
     * Returns a BUGS-language expression representing the function call.
     * The default behaviour for a function named "foo" is to return
     * "foo(arg1,arg2)". Functions that are represented as prefix or infix
     * operators need to override this function.
     *
     * @param par Vector of parameter names for the function
     */
    virtual std::string deparse(std::vector<std::string> const &par) const;
    /**
     * Checks that the required parameters are fixed. Most functions
     * do not require any fixed parameters, so the default
     * implementation returns true. However, for some ArrayFunctions
     * the dimension of the return value depends on the values of some of
     * the parameters. These parameters must be fixed so that the
     * dimension can be determined at compile time. Hence the default
     * function must be overridden for any function with a return
     * value whose dimension cannot be determined purely from the 
     * dimensions of its arguments.
     *
     * @param mask Boolean vector indicating which parameters have
     * fixed values
     */
    virtual bool checkParameterFixed(std::vector<bool> const &mask) const;
    /**
     * Returns true if the function is differentiable with respect to
     * its ith argument. The default implementation returns false.
     * Differentiable functions must override this member function.
     *
     * If this function returns true then it guarantees that the
     * functions ScalarFunction#gradient, VectorFunction#gradient, and
     * ArrayFunction#gradient will correctly calculate the gradient of
     * the function.
     *
     * False negatives are allowed, i.e. the function may return false
     * even if it is differentiable. This means that the corresponding
     * gradient function is not correctly implemented.
     *
     * @param i Index of the parameter to take the derivative
     * (starting from zero).
     */
    virtual bool isDifferentiable(unsigned long i) const;
};

/**
 * Checks that a vector of parameters of length npar is consistent
 * with the function.
 */
inline bool checkNPar(Function const *func, unsigned long npar)
{
    return (func->npar() == 0 && npar > 0) || func->npar() == npar;
}

} /* namespace jags */

#endif /* FUNCTION_H_ */
