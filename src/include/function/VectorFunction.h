#ifndef VECTOR_FUNCTION_H_
#define VECTOR_FUNCTION_H_

#include <function/Function.h>

namespace jags {

/**
 * @short Vector-valued function with vector arguments
 *
 * VectorFunction represents vector-valued functions whose parameters are
 * also vectors.
 */
class VectorFunction : public Function
{
public:
    VectorFunction(std::string const &name, unsigned long npar);
    /**
     * Evaluates the function and writes the result into the given
     * array
     *
     * @param value Array of doubles which will contain the result on exit.
     * @param args Vector of arguments
     * @param lengths Vector of argument lengths: the length of the array
     * of doubles pointed to by args[i] is lengths[i]. 
     */
    virtual void evaluate(double *value, 
			  std::vector <double const *> const &args,
			  std::vector <unsigned long> const &lengths) const = 0;
    /**
     * Calculates the length of the return value based on the arguments.
     *
     * @param arglengths Vector of argument lengths. This must return
     * true when passed to checkParameterLength
     *
     * @param argvalues Vector of pointers to argument values. 
     *
     */
    virtual unsigned long 
	length(std::vector <unsigned long> const &arglengths,
	       std::vector <double const *> const &argvalues) const = 0;
    /**
     * Checks that the lengths of all the arguments are consistent. 
     * The default implementation returns true if all the arguments
     * have non-zero length.
     *
     * This should be overridden by any sub-class that either has
     * additional restrictions on the lengths of the arguments, or may
     * allow zero-length arguments.
     */
    virtual bool 
	checkParameterLength(std::vector<unsigned long> const &args) const;
    /**
     * Checks that the values of the arguments are in the domain of
     * the function. The default implementation always returns true
     * (i.e. assumes there are no restrictions on parameters).
     */
    virtual bool 
	checkParameterValue(std::vector<double const *> const &args,
			    std::vector<unsigned long> const &lengths) const;  
    /**
     * FIXME
     */
    virtual void gradient(double *grad, std::vector<double const *> const &args,
			  std::vector<unsigned long> const &lengths,
			  unsigned long i) const;
};

} /* namespace jags */

#endif /* VECTOR_FUNCTION_H_ */
