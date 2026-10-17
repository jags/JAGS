#ifndef TEST_LIB_H_
#define TEST_LIB_H_

#include <util/logical.h>

#include <vector>

//A mix-in class for test fixtures that provides some useful constants
class JAGSFixture
{
  public:
    JAGSFixture();

  protected:
    //tolerance for equality tests
    double tol; 
    
    //Logical vectors of length 1
    std::vector<bool> F; //False
    std::vector<bool> T; //True
    
    //Logical vectors of length 2
    std::vector<bool> FF; //False, False
    std::vector<bool> FT; 
    std::vector<bool> TF; 
    std::vector<bool> TT; 
    
    //Logical vectors of length 3
    std::vector<bool> FFF; //False, False, False
    std::vector<bool> FFT; 
    std::vector<bool> FTF; 
    std::vector<bool> FTT; 
    std::vector<bool> TFF; 
    std::vector<bool> TFT; 
    std::vector<bool> TTF; 
    std::vector<bool> TTT; 
};

/* Similar to RangeIterator in the JAGS library, this iterates over all
   values of a boolean vector */
class BoolIterator : public std::vector<bool>
{
public:
    /* Construct a boolean iterator of length n with all elements
     * initialized to false */
    BoolIterator(unsigned long n);
    /* Go to the next value */
    void next();
    /* True if we have reached the last value with all elements equal
     * to true */
    bool atEnd;
};

inline bool always(std::vector<bool> const &mask) { return true; }
inline bool never(std::vector<bool> const &mask) { return false; }
inline bool all(std::vector<bool> const &mask) { return jags::allTrue(mask); }
inline bool any(std::vector<bool> const &mask) { return jags::anyTrue(mask); }


/*
  Functions and Distributions in JAGS are set up to take vectors of
  pointers as arguments, with additional arguments for lengths or
  dimensions of the arguments when necessary. This is not a convenient
  format for creating test functions so we use these wrappers to
  convert test function arguments to the correct format.
*/

/* Functions to convert arguments to an STL vector */
   
//Convert a double to a vector of length 1
inline std::vector<double> mkVec(double const &x)
{
    return std::vector<double>(1, x);
}

//Convert a static array of length N to a vector of length N
template<std::size_t N>
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

/* Argument type and return type for array functions and distributions */

typedef std::pair<std::vector<double>, std::vector<unsigned long>> array_value;

/* Functions to convert arguments to array_value
   
 */

/* Convert scalar to array_value */
inline array_value mkArray(double const &x)
{
    return array_value(std::vector<double>(1, x), std::vector<unsigned long>(1, 1UL));
}

/* Convert static array of length N to array_value */
template<std::size_t N>
array_value mkArray(double const (&x)[N])
{
    std::vector<double> y(N);
    copy(x, x + N, y.begin());
    return array_value(y, std::vector<std::size_t>(1, N));
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
    return ans;
}

inline std::vector<double const *> getValues(std::vector<array_value> const &args)
{
    std::vector<double const *> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].first.data();
    }
    return ans;
}

inline std::vector<unsigned long> getLengths(std::vector<std::vector<double>> const &args)
{
    std::vector<unsigned long> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].size();
    }
    return ans;
}

inline std::vector<std::vector<unsigned long>> getDimensions(std::vector<array_value> const &args)
{
    std::vector<std::vector<unsigned long>> ans(args.size());
    for (unsigned long i = 0; i < args.size(); ++i) {
	ans[i] = args[i].second;
    }
    return ans;
}

/* Similar to the R all.equal function we use all_equal to tests approximate
   equality of two objects */

//Test approximate equality of two vectors
bool all_equal(std::vector<double> const &u, std::vector<double> const &v,
	       double tol);

//Test approximate equality of two array_values
bool all_equal(array_value const &A, array_value const &B, double tol);

#endif /* TEST_LIB_H */
