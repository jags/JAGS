#ifndef FUNC_TEST_H_
#define FUNC_TEST_H_

#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "ArrayFunction.h"
#include "LinkFunction.h"

#include <utility>

/*
  Function in JAGS are set up to take vectors of pointers as
  arguments, with additional arguments for lengths or dimensions of
  the arguments when necessary. This is not a good choice for the
  testing framework.
  
  These are wrappers to the members of Function and its subclasses
  that take care of setting up the arguments in the right way. We
  use these to simplify the testing framework.

  The most important functions are "eval" and "veval" for functions
  returning a scalar and a vector, respectively. 

  The eval function can be used to test a ScalarFunction taking from 1
  to 3 arguments.  It takes doubles as arguments and returns a double.

  The veval function can be used to test a VectorFunction taking from
  1 to 4 arguments. Through the use of templates, veval accepts as
  arguments any of the following: an STL vector of doubles, a double,
  or a static array of doubles.

  The eval function is also overloaded to allow evaluation of a
  VectorFunction returning a scalar value.
*/

/* Testing functions valid for all functions */

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

/* Tests for scalar functions */

//Check that the limits of a scalar function are valid
void checkLimits(jags::ScalarFunction const *f, double lower, double upper);

//Evaluate a link function
double eval(jags::LinkFunction const *f, double x);

//Evaluate a scalar function taking a single argument.
double eval(jags::ScalarFunction const *f, double x);
bool checkargs(jags::ScalarFunction const *f, double x);

//Evaluate a scalar function taking two arguments
double eval(jags::ScalarFunction const *f, double x, double y);
bool checkargs(jags::ScalarFunction const *f, double x, double y);

//Evaluate a scalar function taking three arguments
double eval(jags::ScalarFunction const *f, double x, double y, double z);
bool checkargs(jags::ScalarFunction const *f, double x, double y, double z);

//Evaluate gradient for a scalar function taking a single argument
double gradient(jags::ScalarFunction const *f, const double x);
double numgradient(jags::ScalarFunction const *f, const double x, double delta);

//Evaluate gradient for a scalar function taking two arguments
double gradient(jags::ScalarFunction const *f, const double x, double y,
		unsigned long i);
double numgradient(jags::ScalarFunction const *f, const double x, double y,
		   unsigned long i, double delta);

//Evaluate gradient for a scalar function taking three arguments
double gradient(jags::ScalarFunction const *f, const double x, double y,
		double z, unsigned long i);
double numgradient(jags::ScalarFunction const *f, const double x, double y,
		   double z, unsigned long i, double delta);

//Evaluate gradient for a link function
double gradient(jags::LinkFunction const *f, double x);
double numgradient(jags::LinkFunction const *f, double x, double delta);

/* Tests for vector functions */

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

//An apparently trivial conversion function that allows us to mix STL
//vectors with scalars and static arrays as arguments
inline std::vector<double> mkVec(std::vector<double> const &x)
{
    return x;
}

//Safely evaluate a vector function taking a single argument
std::vector<double> veval(jags::VectorFunction const *f, 
			  std::vector<double> const &x);

//Check for valid arguments
bool checkargs(jags::VectorFunction const *f, std::vector<double> const &x);


std::vector<double> vgradient(jags::VectorFunction const *f,
			      std::vector<double> const &x, unsigned long i);

std::vector<double> vnumgradient(jags::VectorFunction const *f,
				 std::vector<double> const &x,
				 unsigned long i, double delta);

//Templated version that allows you to pass any argument that can be
//coerced to an STL vector via mkVec
template<typename T>
std::vector<double> veval(jags::VectorFunction const *f, T const &x)
{
    return veval(f, mkVec(x));
}

template<typename T>
bool checkargs(jags::VectorFunction const *f, T const &x)
{
    return checkargs(f, mkVec(x));
}

template<typename T>
std::vector<double> vgradient(jags::VectorFunction const *f, T const &x,
			      unsigned long i)
{
    return vgradient(f, mkVec(x), i);
}

template<typename T>
std::vector<double> vnumgradient(jags::VectorFunction const *f, T const &x,
				 unsigned long i, double delta)
{
    return vnumgradient(f, mkVec(x), i, delta);
}

//Safely evaluate a vector function taking two arguments
std::vector<double> veval(jags::VectorFunction const *f, 
			  std::vector<double> const &x, 
			  std::vector<double> const &y);

bool checkargs(jags::VectorFunction const *f, 
	       std::vector<double> const &x, 
	       std::vector<double> const &y);

std::vector<double> vgradient(jags::VectorFunction const *f,
			      std::vector<double> const &x,
			      std::vector<double> const &y,
			      unsigned long i);

std::vector<double> vnumgradient(jags::VectorFunction const *f,
				 std::vector<double> const &x,
				 std::vector<double> const &y,
				 unsigned long i, double delta);

//Templated version
template<typename T, typename U>
std::vector<double> veval(jags::VectorFunction const *f,
			  T const &x, U const &y)
{
    return veval(f, mkVec(x), mkVec(y));
}

template<typename T, typename U>
bool checkargs(jags::VectorFunction const *f,
	       T const &x, U const &y)
{
    return checkargs(f, mkVec(x), mkVec(y));
}

template<typename T, typename U>
std::vector<double> vgradient(jags::VectorFunction const *f,
			      T const &x, U const &y,
			      unsigned long i)
{
    return vgradient(f, mkVec(x), mkVec(y), i);
}

template<typename T, typename U>
std::vector<double> vnumgradient(jags::VectorFunction const *f,
				 T const &x, U const &y,
				 unsigned long i, double delta)
{
    return vnumgradient(f, mkVec(x), mkVec(y), i, delta);
}

//Three arguments
std::vector<double> veval(jags::VectorFunction const *f, 
			  std::vector<double> const &x, 
			  std::vector<double> const &y,
			  std::vector<double> const &z);

bool checkargs(jags::VectorFunction const *f, 
	       std::vector<double> const &x, 
	       std::vector<double> const &y,
	       std::vector<double> const &z);

std::vector<double> vgradient(jags::VectorFunction const *f,
			      std::vector<double> const &x,
			      std::vector<double> const &y,
			      std::vector<double> const &z,
			      unsigned long i);

std::vector<double> vnumgradient(jags::VectorFunction const *f,
				 std::vector<double> const &x,
				 std::vector<double> const &y,
				 std::vector<double> const &z,
				 unsigned long i, double delta);

//Three arguments, template
template<typename T, typename U, typename V>
std::vector<double> veval(jags::VectorFunction const *f,
			  T const &x, U const &y, V const &z)
{
    return veval(f, mkVec(x), mkVec(y), mkVec(z));
}


template<typename T, typename U, typename V>
bool checkargs(jags::VectorFunction const *f,
	       T const &x, U const &y, V const &z)
{
    return checkargs(f, mkVec(x), mkVec(y), mkVec(z));
}

template<typename T, typename U, typename V>
std::vector<double> vgradient(jags::VectorFunction const *f,
			      T const &x, U const &y, V const &z,
			      unsigned long i)
{
    return vgradient(f, mkVec(x), mkVec(y), mkVec(z), i);
}

template<typename T, typename U, typename V>
std::vector<double> vnumgradient(jags::VectorFunction const *f,
				 T const &x, U const &y, V const &z,
				 unsigned long i, double delta)
{
    return vnumgradient(f, mkVec(x), mkVec(y), mkVec(z), i, delta);
}

//Four arguments
std::vector<double> veval(jags::VectorFunction const *f, 
			  std::vector<double> const &x, 
			  std::vector<double> const &y,
			  std::vector<double> const &z,
			  std::vector<double> const &w);

bool checkargs(jags::VectorFunction const *f, 
	       std::vector<double> const &x, 
	       std::vector<double> const &y,
	       std::vector<double> const &z,
	       std::vector<double> const &w);

//Four arguments, template
template<typename T1, typename T2, typename T3, typename T4>
std::vector<double> veval(jags::VectorFunction const *f, T1 const &x1,
			  T2 const &x2, T3 const &x3, T4 const &x4)
{
    return veval(f, mkVec(x1), mkVec(x2), mkVec(x3), mkVec(x4));
}

template<typename T1, typename T2, typename T3, typename T4>
bool checkargs(jags::VectorFunction const *f,
	       T1 const &x1, T2 const &x2, T3 const &x3, T4 const &x4)
{
    return checkargs(f, mkVec(x1), mkVec(x2), mkVec(x3), mkVec(x4));
}


/* Vector functions returning a scalar */

//Single argument
double eval(jags::VectorFunction const *f, std::vector<double> const &x);

//Template version
template <typename T>
double eval(jags::VectorFunction const *f, T const &x)
{
    return eval(f, mkVec(x));
}

//Two arguments
double eval(jags::VectorFunction const *f, std::vector<double> const &x, 
	    std::vector<double> const &y);

//Template version
template<typename T, typename U>
double eval(jags::VectorFunction const *f, T const &x, U const &y)
{
    return eval(f, mkVec(x), mkVec(y));
}

//Three arguments
double eval(jags::VectorFunction const *f, std::vector<double> const &x, 
	    std::vector<double> const &y, std::vector<double> const &z);

//Template version
template<typename T, typename U, typename V>
double eval(jags::VectorFunction const *f, T const &x, U const &y, V const &z)
{
    return eval(f, mkVec(x), mkVec(y), mkVec(z));
}

/* Argument type and return type for aeval */

typedef std::pair<std::vector<double>, std::vector<unsigned long>> array_value;

//Test approximate equality of two array_values
bool all_equal(array_value const &A, array_value const &B, double tol);

/* Array function taking a single argument */

//Safely evaluate an array function taking a single argument
array_value aeval(jags::ArrayFunction const *f, array_value const &x);

bool checkargs(jags::ArrayFunction const *f, array_value const &x);

/* Array function taking two arguments */

array_value aeval(jags::ArrayFunction const *f, array_value const &x,
		  array_value const &y);

bool checkargs(jags::ArrayFunction const *f, array_value const &x,
	       array_value const &y);


// Array function taking 3 arguments

array_value aeval(jags::ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z);

bool checkargs(jags::ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z);

// Array function taking 4 arguments

array_value aeval(jags::ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z,
		  array_value const &u);

bool checkargs(jags::ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z,
	       array_value const &u);

// Array function taking 5 arguments

array_value aeval(jags::ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z,
		  array_value const &u, array_value const &v);

bool checkargs(jags::ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z,
	       array_value const &u, array_value const &v);

/* Array functions returning a scalar */

double eval(jags::ArrayFunction const *f, array_value const &x);

double eval(jags::ArrayFunction const *f, array_value const &x,
	    array_value const &y);

double eval(jags::ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z);

double eval(jags::ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z,
	    array_value const &u);

double eval(jags::ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z,
	    array_value const &u, array_value const &v);

#endif /* FUNC_TEST_H_ */
