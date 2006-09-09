#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <string>
#include <vector>

class SArray;
class Index;
class Range;

/**
 * An abstract base class for functions.
 *
 * Function objects contain only constant data members and all member
 * functions are constant. Hence only one object needs to be
 * instantiated for each subclass.
 *
 * @see FuncTab
 * @short Function
 */
class Function
{
  const std::string _name;
  const unsigned int _npar;
 public:
  /**
   * Constructor.
   * @param name name of the function as used in the BUGS language
   * @param npar number of parameters. If npar == 0 then the function
   * takes a variable number of parameters.
   */
  Function (std::string const &name, unsigned int npar);
  virtual ~Function ();
  /**
   * Returns the  BUGS language name of the function
   */
  std::string const &name () const;
  /**
   * Evaluates the function. 
   *
   * @param value SArray which is modified to contain the result of
   * the evaluation 
   * @param args Vector of arguments. This is assumed to satisfy
   * the requirements of length, dimension, and compatibility with
   * the dimension of value.
   */
  virtual void
    evaluate (SArray & value,
	      std::vector <SArray const *> const &args) const = 0;
  /**
   * Checks that a vector of parameters of length npar is consistent
   * with the function.
   */
  bool checkParameterLength(unsigned int npar) const;
  /**
   * Checks whether dimensions of the function parameters are correct.  
   *
   * @param dims Vector of length npar denoting the dimensions of
   * the parameters, with any redundant dimensions dropped.
   */
  virtual bool
     checkParameterDim (std::vector <Index> const &dims) const = 0;
  /**
   * Checks whether the parameter values lie in the domain of the
   * function.  
   *
   * Since many functions do not need to check their parameter
   * values, a default method is provided which always returns
   * true.
   */
  virtual bool
    checkParameterValue (std::vector <SArray const *> const &args) const;
  /**
   * Calculates what the dimension of the return value should be,
   * based on the dimensions of the arguments. The dimension of the
   * return value cannot depend on the value of any of the arguments,
   * only their dimensions.
   *
   * @param dims Vector of Indices denoting the dimensions of the 
   * parameters. This vector must return true when passed to
   * checkParameterDim
   */
  virtual Index dim (std::vector <Index> const &dims) const = 0;
  /**
   * Returns true if the function returns integer values. The default
   * implementation returns false, so any function that can return
   * integer values will need to overload this.
   *
   * @param mask Vector indicating whether parameters are discrete
   * or not. 
   *
   * @see SArray##isDiscreteValued
   */
  virtual bool isDiscreteValued(std::vector<bool> const &mask) const;
  /**
   * Tests whether the value of the function is linear.  False
   * negative responses are permitted: i.e. the value false may be
   * returned when the value is, in fact, a linear function of the
   * given parameters, but false positives are not allowed. Since most
   * functions of interest are non-linear, the default method returns
   * false, so linear functions must overload this.
   * 
   * @param mask boolean vector of length equal to the number of
   * parameters.  The mask indicates a subset of parameters (those with
   * value true) to be tested for linearity. The function may be 
   * non-linear for parameters with mask value false.
   *
   * @param isfixed boolean vector. This may be empty, in which case
   * it is ignored, and the function tests for linearity of the
   * parameters in mask. A non-empty vector must have length equal to
   * mask, and denotes the parameters whose values are fixed. In this
   * case the test is more restrictive, and the function returns true
   * only if the parameters of the coefficients are fixed, i.e
   * the function is of the form f(x) = A + B * x where B is fixed
   * (but A may vary).
   *
   * @seealso DeterministicNode##isLinear
   */
  virtual bool isLinear(std::vector<bool> const &mask,
			std::vector<bool> const &isfixed) const;
  /**
   * Tests whether the value of the function is a scale transformation
   * of the given parameter.  False negative responses are permitted.
   * Since most functions of interest are not scale transormations,
   * the default method returns false.
   * 
   * @param index The index of the parameter to be tested. The function
   * may be a non-scale transformation of the other pararameters (e.g.
   * the function x/y is a scale transformation of x but not of y).
   *
   * @param isfixed boolean vector. This may be empty, in which case
   * it is ignored.  A non-empty vector must have length equal to
   * mask, and denotes the parameters whose values are fixed. In this
   * case the test is more restrictive, and the function returns true
   * only if the scale transformation is fixed, i.e. the function is
   * of the form f(x) = B * x where B is fixed.
   *
   * @seealso DeterministicNode##isScale
   */
  virtual bool isScale(unsigned int index,
		       std::vector<bool> const &isfixed) const;
};

#endif /* FUNCTION_H_ */
