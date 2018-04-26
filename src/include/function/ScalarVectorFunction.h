#ifndef SCALAR_VECTOR_FUNCTION_H_
#define SCALAR_VECTOR_FUNCTION_H_

#include <function/VectorFunction.h>

namespace jags {

    /**
     * @short Scalar-valued function with vector arguments
     *
     * ScalarVectorFunction represents scalar-valued functions whose
     * parameters are vectors.
     */
    class ScalarVectorFunction : public VectorFunction
    {
    public:
	ScalarVectorFunction(std::string const &name, unsigned long npar);
	/**
	 * Evaluates the function and returns the scalar result
	 * @param args Vector of arguments
	 * @param lengths Vector of argument lengths: the length of the array
	 * of doubles pointed to by args[i] is lengths[i]. 
	 */
	virtual double 
	    scalarEval (std::vector <double const *> const &args,
			std::vector <unsigned long> const &lengths) const = 0;
	/**
	 * Implements VectorFunction#evaluate in terms of 
	 * ScalarVectorFunction#eval.
	 */
	void evaluate(double *value, 
		      std::vector <double const *> const &args,
		      std::vector <unsigned long> const &lengths) const;
	/**
	 * Implements VectorFunction#length, returning the value 1
	 */
	unsigned long 
	    length(std::vector <unsigned long> const &arglengths,
		   std::vector <double const *> const &argvalues) const;
    };

} /* namespace jags */

#endif /* SCALAR_VECTOR_FUNCTION_H_ */
