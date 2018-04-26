#ifndef ARRAY_FUNCTION_H_
#define ARRAY_FUNCTION_H_

#include <function/Function.h>

#include <string>
#include <vector>

namespace jags {

/**
 * @short Matrix- or array-valued function
 *
 * Array-valued functions are the most general class of function. The
 * arguments of an array-valued function, and the value may be a
 * scalar, vector, or array.
 *
 */
class ArrayFunction : public Function
{
public:
    /**
     * Constructor.
     */
    ArrayFunction(std::string const &name, unsigned long npar);
    /**
     * Evaluates the function. 
     *
     * @param value array of doubles which contains the result of
     * the evaluation on exit 
     * @param args Vector of arguments. 
     * @param dims Respective dimensions of each element of args.
     */
    virtual void
	evaluate(double *value, 
		 std::vector<double const *> const &args,
		 std::vector<std::vector<unsigned long> > const &dims) const = 0;
    /**
     * Checks whether dimensions of the function parameters are correct.  
     *
     * @param dims Vector of length npar denoting the dimensions of
     * the parameters, with any redundant dimensions dropped.
     */
    virtual bool
	checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) 
	const = 0;
    /**
     * Checks whether the parameter values lie in the domain of the
     * function. The default implementation returns true.
     */
    virtual bool
	checkParameterValue(std::vector<double const *> const &args,
                            std::vector<std::vector<unsigned long> > const &dims) 
        const;
    /**
     * Calculates what the dimension of the return value should be,
     * based on the arguments. 
     *
     * @param dims Vector of Indices denoting the dimensions of the
     * parameters. This vector must return true when passed to
     * checkParameterDim.
     *
     * @param values Vector of pointers to parameter values.
     */
    virtual std::vector<unsigned long> 
	dim(std::vector <std::vector<unsigned long> > const &dims,
	    std::vector <double const *> const &values) const = 0;
};

} /* namespace jags */

#endif /* ARRAY_FUNCTION_H_ */
