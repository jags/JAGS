#ifndef FUNC_TEST_H_
#define FUNC_TEST_H_

#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "ArrayFunction.h"
#include "LinkFunction.h"

#include <utility>

//For transform. FIXME: Move this intot testfuncc
#include <algorithm>

/*
  Function in JAGS are set up to take vectors of pointers as
  arguments, with additional arguments for lengths or dimensions of
  the arguments when necessary. 
  
  These are wrappers to the members of Function and its subclasses
  that simplify the testing framework. 

  Variadic function templates are used to pass the arguments to the
  testing functions in the testing library. Arguments may be scalars,
  STL vectors, static arrays or, for array functions, a special struct
  called array_value. Conversion functions are used to convert the
  arguments to the appropriate format. Forward references are used
  to pass the arguments on to the conversion functions, so that we may
  supply both lvalues and rvalues.
*/

/* Functions to convert arguments to an STL vector */
   
//Convert a double to a vector of length 1
inline std::vector<double> mkVec(double const &x)
{
    return std::vector<double>(1, x);
}

//Convert a static array of length N to a vector of length N
template<size_t N>
std::vector<double> mkVec(double const (&x)[N])
{
    std::vector<double> y(N);
    copy(x, x + N, y.begin());
    return y;
}

//Return a copy of a vector
inline std::vector<double> mkVec(std::vector<double> const &x)
{
    return x;
}

/* Argument type and return type for array functions */

typedef std::pair<std::vector<double>, std::vector<unsigned long>> array_value;

/* Functions to convert arguments to array_value
   
 */

/* Convert scalar to array_value */
inline array_value mkArray(double const &x)
{
    return array_value(std::vector<double>(1, x), std::vector<unsigned long>(1, 1UL));
}

/* Convert static array of length N to array_value */
template<size_t N>
array_value mkArray(double const (&x)[N])
{
    std::vector<double> y(N);
    copy(x, x + N, y.begin());
    return array_value(y, std::vector<size_t>(1, N));
}

//Convert vector to array_value
inline array_value mkArray(std::vector<double> const &x)
{
    return array_value(x, std::vector<unsigned long>(1, x.size()));
}

//Return a copy of an array_value
inline array_value mkArray(array_value const &x)
{
    return x;
}

/* After converting the arguments to the appropriate class, we need to extract values,
   lengths, dimensions, ... in a format suitable for passing on to the JAGS library */

inline std::vector<double const *> getValues(std::vector<std::vector<double>> const &args)
{
    std::vector<double const *> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].data();
    }
    //FIXME GET THIS WORKING
    //std::transform(args.begin(), args.end(), ans.begin(), getData);
    return ans;
}

inline std::vector<double const *> getValues(std::vector<array_value> const &args)
{
    std::vector<double const *> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].first.data();
    }
    //std::transform(args.begin(), args.end(), ans.begin(), getData);
    return ans;
}

inline std::vector<unsigned long> getLengths(std::vector<std::vector<double>> const &args)
{
    std::vector<unsigned long> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].size();
    }
    //std::transform(args.begin(), args.end(), ans.begin(), getLen);
    return ans;
}

inline std::vector<std::vector<unsigned long>> getDimensions(std::vector<array_value> const &args)
{
    std::vector<std::vector<unsigned long>> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].second;
    }
    //std::transform(args.begin(), args.end(), ans.begin(), getLen);
    return ans;
}


/* Check function arguments */

/* Scalar functions */

/* Check that arguments of a scalar function are valid */
bool checkArgs(jags::ScalarFunction const *f,
	       std::vector<double const *> const &args);

/* Variadic template version of checkArgs */
template<typename... Args>
bool checkargs(jags::ScalarFunction const *f, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return checkArgs(f, getValues(vargs));
}

/* Check that the limits of a scalar function, i.e. that lower and
   upper are valid arguments, but a value below lower or above upper
   is invalid */
void checkLimits(jags::ScalarFunction const *f, double lower, double upper);

/* Vector functions */

/* Check that argumetns of a vector function are valid */
bool checkArgs(jags::VectorFunction const *f,
	       std::vector<double const *> const &args,
	       std::vector<unsigned long> const &arglen);

/* Variadic template version of checkArgs for vector functions */
template<typename... Args>
bool checkargs(jags::VectorFunction const *f, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return checkArgs(f, getValues(vargs), getLengths(vargs));
}

// Array functions

bool checkArgs(jags::ArrayFunction const *f,
	       std::vector<double const *> const &args,
	       std::vector<std::vector<unsigned long>> const &argdims);

template<typename... Args>
bool checkargs(jags::ArrayFunction const *f, Args&&... args)
{
    std::vector<array_value> aargs{mkArray(std::forward<Args>(args))...};
    return checkArgs(f, getValues(aargs), getDimensions(aargs));
}

/* Safely evaluate functions after applying argument checks */

//Link functions

double eval(jags::LinkFunction const *f, double x);

// Scalar functions

double Eval(jags::ScalarFunction const *f,
	    std::vector<double const *> const &args);

template<typename... Args>
double eval(jags::ScalarFunction const *f, Args&&... args)
{
    //std::vector<std::vector<double>> vargs{mkVec(args)...};
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return Eval(f, getValues(vargs));
}

// Vector functions

std::vector<double> VEval(jags::VectorFunction const *f,
			  std::vector<double const *> const &args,
			  std::vector<unsigned long> const &arglen);

template<typename... Args>
std::vector<double> veval(jags::VectorFunction const *f, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return VEval(f, getValues(vargs), getLengths(vargs));
}

// Array functions

array_value AEval(jags::ArrayFunction const *f,
		  std::vector<double const *> const &args,
		  std::vector<std::vector<unsigned long>> const &argdims);

template<typename... Args>
array_value aeval(jags::ArrayFunction const *f, Args&&... args)
{
    std::vector<array_value> aargs{mkArray(std::forward<Args>(args))...};
    return AEval(f, getValues(aargs), getDimensions(aargs));
}

// Convenience wrappers for vector and array functions returning a scalar

double Eval(jags::VectorFunction const *f,
	    std::vector<double const *> const &args,
	    std::vector<unsigned long> const &arglen);

template<typename... Args>
double eval(jags::VectorFunction const *f, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return Eval(f, getValues(vargs), getLengths(vargs));
}

double Eval(jags::ArrayFunction const *f,
	    std::vector<double const *> const &args,
	    std::vector<std::vector<unsigned long>> const &argdims);

template<typename... Args>
double eval(jags::ArrayFunction const *f, Args&&... args)
{
    std::vector<array_value> aargs{mkArray(std::forward<Args>(args))...};
    return Eval(f, getValues(aargs), getDimensions(aargs));
}

/* Testing functions valid for all function classes */

//Check all possible values of mask using a predicate (see below)
bool isdiscrete(jags::Function const *f, unsigned long npar,
		bool (*predicate) (std::vector<bool> const &));

//suitable predicates for isdiscrete
bool always(std::vector<bool> const &); //returns true
bool never(std::vector<bool> const &); //returns false
bool all(std::vector<bool> const &); //returns true if all arguments are true
bool any(std::vector<bool> const &); //returns true if any argyments are true

//Returns true if f is never an additive function
bool neveradditive(jags::Function const *f, unsigned long npar);
//Returns true if f is never a scale function
bool neverscale(jags::Function const *f, unsigned long npar);
//Returns true if f is never a linear function
bool neverlinear(jags::Function const *f, unsigned long npar);
//Returns true if f is never a power function
bool neverpow(jags::Function const *f, unsigned long npar);
//Returns true if f is never an additive, linear, scale, or power function
bool neverclosed(jags::Function const *f, unsigned long npar);

/*  Evaluate gradients with parameter checks */

/* Link functions */

double gradient(jags::LinkFunction const *f, double x);
double numgradient(jags::LinkFunction const *f, double x, double delta);

/* Scalar functions */

double Gradient(jags::ScalarFunction const *f,
		std::vector<double const *> const &args,
		unsigned long i);

double NumGradient(jags::ScalarFunction const *f,
		   std::vector<double const*> const &args,
		   unsigned long i, double delta);

template<typename... Args>
double gradient(jags::ScalarFunction const *f, unsigned long i, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return Gradient(f, getValues(vargs), i);
}

template<typename... Args>
double numgradient(jags::ScalarFunction const *f, unsigned long i, double delta, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return NumGradient(f, getValues(vargs), i, delta);
}

/* Convenience wrappers for scalar functions with one argument */

inline double gradient(jags::ScalarFunction const *f, double x)
{
    return gradient(f, 0UL, x);
}

inline double numgradient(jags::ScalarFunction const *f, double x, double delta)
{
    return numgradient(f, 0UL, delta, x);
}

/* Vector functions */

std::vector<double> VGradient(jags::VectorFunction const *f,
			      std::vector<double const *> const &args,
			      std::vector<unsigned long> const &arglens,
			      unsigned long i);

template<typename... Args>
std::vector<double> vgradient(jags::VectorFunction const *f, unsigned long i, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return VGradient(f, getValues(vargs), getLengths(vargs), i);
}


std::vector<double> VNumGradient(jags::VectorFunction const *f,
				 std::vector<double const*> const &args,
				 std::vector<unsigned long> const &arglen,
				 unsigned long i, double delta);

template<typename... Args>
std::vector<double> vnumgradient(jags::VectorFunction const *f, unsigned long i, double delta, Args&&... args)
{
    std::vector<std::vector<double>> vargs{mkVec(std::forward<Args>(args))...};
    return VNumGradient(f, getValues(vargs), getLengths(vargs), i, delta);
}

// Array functions

std::vector<double> VGradient(jags::ArrayFunction const *f,
			      std::vector<double const *> const &args,
			      std::vector<std::vector<unsigned long>> const &arglens,
			      unsigned long i);

template<typename... Args>
std::vector<double> vgradient(jags::ArrayFunction const *f, unsigned long i, Args&&... args)
{
    std::vector<array_value> aargs{mkArray(std::forward<Args>(args))...};
    return VGradient(f, getValues(aargs), getDimensions(aargs), i);
}


std::vector<double> VNumGradient(jags::ArrayFunction const *f,
				 std::vector<double const*> const &args,
				 std::vector<std::vector<unsigned long>> const &arglen,
				 unsigned long i, double delta);

template<typename... Args>
std::vector<double> vnumgradient(jags::ArrayFunction const *f, unsigned long i, double delta, Args&&... args)
{
    std::vector<array_value> aargs{mkArray(std::forward<Args>(args))...};
    return VNumGradient(f, getValues(aargs), getDimensions(aargs), i, delta);
}

//Test approximate equality of two vectors
bool all_equal(std::vector<double> const &u, std::vector<double> const &v,
	       double tol);

//Test approximate equality of two array_values
bool all_equal(array_value const &A, array_value const &B, double tol);

#endif /* FUNC_TEST_H_ */
