#include <function/testfun.h>
#include <util/logical.h>
#include <util/dim.h>

#include <cppunit/extensions/HelperMacros.h>

using jags::ScalarFunction;
using jags::VectorFunction;
using jags::ArrayFunction;
using jags::LinkFunction;
using jags::Function;
using jags::anyTrue;
using jags::allTrue;
using jags::product;

#include <climits>
#include <cmath>
#include <algorithm>

using std::vector;
using std::string;
using std::copy;
using std::floor;
using std::isfinite;
using std::isnan;
using std::pair;

/* All functions */

bool isdiscrete(Function const *f, bool mask1)
{
    CPPUNIT_ASSERT(checkNPar(f, 1));
    return f->isDiscreteValued(vector<bool>(1, mask1));
}

bool isdiscrete(Function const *f, bool mask1, bool mask2)
{
    CPPUNIT_ASSERT(checkNPar(f, 2));
    vector<bool> arg(2);
    arg[0] = mask1;
    arg[1] = mask2;
    return f->isDiscreteValued(arg);
}

bool isdiscrete(Function const *f, bool mask1, bool mask2, bool mask3)
{
    CPPUNIT_ASSERT(checkNPar(f, 3));
    vector<bool> arg(3);
    arg[0] = mask1;
    arg[1] = mask2;
    arg[2] = mask3;
    return f->isDiscreteValued(arg);
}

class BoolIterator : public std::vector<bool>
{
public:
    bool atEnd;
    BoolIterator(unsigned long n) : vector<bool>(n, false), atEnd(false) {};

    void next() {
	bool bump = true;
	for (unsigned long i = 0; i < size(); ++i) {
	    if (bump) {
		bool x = operator[](i); //current value
		bump = x;
		operator[](i) = !x;
	    }
	    else return;
	}
	if (bump) atEnd=true;
    }
    
};

bool isdiscrete(Function const *f, unsigned long npar,
		bool (*predicate) (vector<bool> const &))
{
    CPPUNIT_ASSERT(checkNPar(f, npar));
    BoolIterator mask(npar);
    
    for(BoolIterator mask(npar); !mask.atEnd; mask.next()) {
	if (f->isDiscreteValued(mask) != predicate(mask)) {
	    return false;
	}
    }
    return true;
}

bool always(vector<bool> const &mask) { return true; }
bool never(vector<bool> const &mask) { return false; }
bool all(vector<bool> const &mask) { return allTrue(mask); }
bool any(vector<bool> const &mask) { return anyTrue(mask); }

bool neveradditive(Function const *f, unsigned long npar)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkNPar(f, npar));

    for(BoolIterator mask(npar); !mask.atEnd; mask.next()) {
	if (f->isAdditive(mask, vector<bool>())) return false;
	for(BoolIterator fixed(npar); !fixed.atEnd; fixed.next()) {
	    if (f->isAdditive(mask, fixed)) return false;
	}
    }
    return true;
}

bool neverlinear(Function const *f, unsigned long npar)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkNPar(f, npar));

    for(BoolIterator mask(npar); !mask.atEnd; mask.next()) {
	if (f->isLinear(mask, vector<bool>())) return false;
	for(BoolIterator fixed(npar); !fixed.atEnd; fixed.next()) {
	    if (f->isLinear(mask, fixed)) return false;
	}
    }
    return true;
}

bool neverscale(Function const *f, unsigned long npar)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkNPar(f, npar));

    for(BoolIterator mask(npar); !mask.atEnd; mask.next()) {
	if (f->isScale(mask, vector<bool>())) return false;
	for(BoolIterator fixed(npar); !fixed.atEnd; fixed.next()) {
	    if (f->isScale(mask, fixed)) return false;
	}
    }
    return true;
}

bool neverpow(Function const *f, unsigned long npar)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkNPar(f, npar));

    for(BoolIterator mask(npar); !mask.atEnd; mask.next()) {
	if (f->isPower(mask, vector<bool>())) return false;
	for(BoolIterator fixed(npar); !fixed.atEnd; fixed.next()) {
	    if (f->isPower(mask, fixed)) return false;
	}
    }
    return true;
}

bool neverclosed(Function const *f, unsigned long npar)
{
    return neverscale(f, npar) && neverlinear(f, npar) &&
	neverpow(f, npar) && neveradditive(f, npar);
}

/* Scalar functions */

static bool checkval(ScalarFunction const *f, double x)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkNPar(f, 1));
    vector<double const *> arg(1, &x);
    return f->checkParameterValue(arg);
}

void checkLimits(ScalarFunction const *f, double lower, double upper)
{
    CPPUNIT_ASSERT(lower < upper);
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkval(f, lower));
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkval(f, upper));

    CPPUNIT_ASSERT_MESSAGE(f->name(), !isnan(eval(f, lower)));
    CPPUNIT_ASSERT_MESSAGE(f->name(), !isnan(eval(f, upper)));

    if (isfinite(upper)) {
	if (upper > 1) {
	    upper *= (1.0 + DBL_EPSILON);
	}
	else if (upper < -1) {
	    upper *= (1.0 - DBL_EPSILON);
	}
	else {
	    upper += DBL_EPSILON;
	}
	CPPUNIT_ASSERT_MESSAGE(f->name(), !checkval(f, upper));
    }
    if (isfinite(lower)) {
	if (lower > 1) {
	    lower *= (1.0 - DBL_EPSILON);
	}
	else if (lower < -1) {
	    lower *= (1.0 + DBL_EPSILON);
	}
	else {
	    lower -= DBL_EPSILON;
	}
	CPPUNIT_ASSERT_MESSAGE(f->name(), !checkval(f, lower));
    }
}

    
static vector<bool> discreteMask(vector<double const *> const &args)
{
    vector<bool> out(args.size(), true);
    
    for (unsigned long i = 0; i < args.size(); ++i) {
	double v = *args[i];
	if (v != floor(v + 0.5)) {
	    out[i] = false;
	}
    }

    return out;
}

static bool checkArgs(ScalarFunction const *f,
		      vector<double const *> const &args)
{
    //Check that arguments are valid
    return checkNPar(f, args.size()) &&
	f->checkParameterDiscrete(discreteMask(args)) &&
	f->checkParameterValue(args);
}

static double checkEval(ScalarFunction const *f,
			vector<double const *> const &args)
{
    //Evaluate scalar function with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args));
    return f->evaluate(args);
}

static double checkGradient(ScalarFunction const *f,
			    vector<double const *> const &args,
			    unsigned long i)
{
    //Evaluate gradient with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));
    return f->gradient(args, i);
}

static double numericGradient(ScalarFunction const *f,
			      vector<double const*> const &args,
			      unsigned long i, double delta)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));

    //Create mutable copy of the arguments
    unsigned long N = args.size();
    vector<double> args0(N);
    vector<double const *> args1(N);
    for (unsigned long j = 0; j < N; ++j) {
	args0[j] = *args[j];
	args1[j] = &args0[j];
    }

    args0[i] = *args[i] - delta;
    double y1 = checkEval(f, args1);
    args0[i] = *args[i] + delta;
    double y2 = checkEval(f, args1);
	
    return (y2 - y1)/(2*delta);
}

/* 
   Evaluate a scalar function that takes a single argument
*/

static vector<double const *> mkArgs(double const *x)
{
    return vector<double const *>(1, x);
}

double eval(ScalarFunction const *f, const double x)
{
    return checkEval(f, mkArgs(&x));
}

bool checkargs(ScalarFunction const *f, const double x)
{
    return checkArgs(f, mkArgs(&x));
}

double gradient(ScalarFunction const *f, const double x)
{
    return checkGradient(f, mkArgs(&x), 0);
}

double numgradient(ScalarFunction const *f, const double x, double delta)
{
    return numericGradient(f, mkArgs(&x), 0, delta);
}

/* Link functions */

double eval(LinkFunction const *f, double x)
{
    return f->inverseLink(x);
}

double gradient(LinkFunction const *f, double x)
{
    return f->grad(x);
}

double numgradient(LinkFunction const *f, double x, double delta)
{
    double y1 = eval(f, x - delta);
    double y2 = eval(f, x + delta);
	
    return (y2 - y1)/(2*delta);
}

/* 
   Evaluate a scalar function that takes two arguments
*/
vector<double const *> mkArgs(double const *x, double const *y)
{
    vector<double const *> args(2);
    args[0] = x;
    args[1] = y;
    return args;
}

double eval(ScalarFunction const *f, double x, double y)
{
    return checkEval(f, mkArgs(&x, &y));
}

bool checkargs(ScalarFunction const *f, double x, double y)
{
    return checkArgs(f, mkArgs(&x, &y));
}

double gradient(ScalarFunction const *f, const double x, double y,
		unsigned long i)
{
    CPPUNIT_ASSERT(i < 2UL);
    return checkGradient(f, mkArgs(&x, &y), i);
}

double numgradient(ScalarFunction const *f, const double x, double y,
		   unsigned long i, double delta)
{
    CPPUNIT_ASSERT(i < 2UL);
    return numericGradient(f, mkArgs(&x, &y), i, delta);
}

/* 
   Evaluate a scalar function that takes three arguments
*/

vector<double const *> mkArgs(double const *x, double const *y, double const *z)
{
    vector<double const *> args(3);
    args[0] = x;
    args[1] = y;
    args[2] = z;
    return args;
}

double eval(ScalarFunction const *f, double x, double y, double z)
{
    return checkEval(f, mkArgs(&x, &y, &z));
}

bool checkargs(ScalarFunction const *f, double x, double y, double z)
{
    return checkArgs(f, mkArgs(&x, &y, &z));
}

double gradient(ScalarFunction const *f, const double x, double y, double z,
		unsigned long i)
{
    CPPUNIT_ASSERT(i < 3UL);
    return checkGradient(f, mkArgs(&x, &y, &z), i);
}

double numgradient(ScalarFunction const *f, const double x, double y, double z,
		   unsigned long i, double delta)
{
    CPPUNIT_ASSERT(i < 3UL);
    return numericGradient(f, mkArgs(&x, &y, &z), i, delta);
}

static vector<bool> discreteMask(vector<double const *> const &args,
				 vector<unsigned long> const &arglen)
{
    vector<bool> out(args.size(), true);
    
    for (unsigned long i = 0; i < args.size(); ++i) {
	double const *v = args[i];
	for (unsigned long j = 0; j < arglen[i]; ++j) {
	    if (v[j] != floor(v[j] + 0.5)) {
		out[i] = false;
		break;
	    }
	}
    }

    return out;
}
			   
static bool checkVArgs(VectorFunction const *f,
		       vector<double const *> const &args,
		       vector<unsigned long> const &arglen)
{
    return args.size() == arglen.size() &&
	checkNPar(f, args.size()) &&
	f->checkParameterLength(arglen) &&
	f->checkParameterDiscrete(discreteMask(args, arglen)) &&
	f->checkParameterValue(args, arglen);
}

static vector<double> checkVEval(VectorFunction const *f,
				 vector<double const *> const &args,
				 vector<unsigned long> const &arglen)
{
    // Evaluate vector function with checks
    CPPUNIT_ASSERT_MESSAGE(string("Valid arguments for ") + f->name(),
			   checkVArgs(f, args, arglen));
    vector<double> ans(f->length(arglen, args));
    f->evaluate(&ans[0], args, arglen);
    return ans;
}

static vector<double> checkVGrad(VectorFunction const *f,
				 vector<double const *> const &args,
				 vector<unsigned long> const &arglen,
				 unsigned long i)
{
    //Evaluate gradient with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkVArgs(f, args, arglen));
    unsigned long n = f->length(arglen, args);
    unsigned long m = arglen[i];
    vector<double> ans(n * m, 0);

    f->gradient(ans.data(), args, arglen, i);
    return ans;
}

static vector<double> numericVGrad(VectorFunction const *f,
				   vector<double const *> const &args,
				   vector<unsigned long> const &arglen,
				   unsigned long i, double delta)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkVArgs(f, args, arglen));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));

    //Create mutable copy of the arguments
    unsigned long N = args.size();
    vector<vector<double>> args0(N);
    vector<double const *> args1(N);
    for  (unsigned long i = 0; i < args.size(); ++i) {
	args0[i] = vector<double>(arglen[i]);
	copy(args[i], args[i] + arglen[i], args0[i].begin());
	args1[i] = args0[i].data();
    }

    //Dimensions of answer matrix
    unsigned long n = f->length(arglen, args);
    unsigned long m = arglen[i];

    vector<double> ans(n * m, 0);
    for (unsigned long j = 0; j < arglen[i]; ++j) {
	args0[i][j] = args[i][j] - delta;
	vector<double> y1 = checkVEval(f, args, arglen);
	args0[i][j] = args[i][j] + delta;
	vector<double> y2 = checkVEval(f, args, arglen);
	args0[i][j] = args[i][j];
	for (unsigned long k = 0; k < n; ++k) {
	    ans[j*m + k] = (y2[k] - y1[k])/(2*delta);
	    //ans[k*n + j] = (y2[k] - y1[k])/(2*delta); ??
	}
    }

    return ans;
}

/* Evaluate a VectorFunction that takes a single argument */
static vector<double const *> mkArgs(vector<double> const &x)
{
    return vector<double const *>(1, &x[0]);
}

static vector<unsigned long> mkLens(vector<double> const &x)
{
    return vector<unsigned long>(1, x.size());
}

vector<double> veval(VectorFunction const *f, vector<double> const &x)
{
    return checkVEval(f, mkArgs(x), mkLens(x));
}

bool checkargs(VectorFunction const *f, vector<double> const &x)
{
    return checkVArgs(f, mkArgs(x), mkLens(x));
}

vector<double> vgradient(VectorFunction const *f, vector<double> const &x,
			unsigned long i)
{
    CPPUNIT_ASSERT(i == 0);
    return checkVGrad(f, mkArgs(x), mkLens(x), i);
}

vector<double> vnumgradient(VectorFunction const *f, vector<double> const &x,
			   unsigned long i, double delta)
{
    CPPUNIT_ASSERT(i == 0);
    return numericVGrad(f, mkArgs(x), mkLens(x), i, delta);
}

/* Evaluate a VectorFunction that takes two arguments */

static vector<double const *>
mkArgs(vector<double> const &x, vector<double> const &y)
{
    vector<double const *> arg(2);
    arg[0] = &x[0];
    arg[1] = &y[0];
    return arg;
}

static vector<unsigned long>
mkLens(vector<double> const &x, vector<double> const &y)
{
    vector<unsigned long> arglen(2);
    arglen[0] = x.size();
    arglen[1] = y.size();
    return arglen;
}

vector<double>
veval(VectorFunction const *f, vector<double> const &x, vector<double> const &y)
{
    return checkVEval(f, mkArgs(x,y), mkLens(x,y));
}

bool checkargs(VectorFunction const *f, vector<double> const &x,
	       vector<double> const &y)
{
    return checkVArgs(f, mkArgs(x,y), mkLens(x,y));
}

vector<double> vgradient(VectorFunction const *f, vector<double> const &x,
			 vector<double> const &y, unsigned long i)
{
    CPPUNIT_ASSERT(i < 2UL);
    return checkVGrad(f, mkArgs(x, y), mkLens(x, y), i);
}

vector<double> vnumgradient(VectorFunction const *f, vector<double> const &x,
			    vector<double> const &y,
			    unsigned long i, double delta)
{
    CPPUNIT_ASSERT(i < 2UL);
    return numericVGrad(f, mkArgs(x, y), mkLens(x, y), i, delta);
}


static vector<double const *> mkArgs(vector<double> const &x,
				     vector<double> const &y,
				     vector<double> const &z)
{
    vector<double const *> arg(3);
    arg[0] = &x[0];
    arg[1] = &y[0];
    arg[2] = &z[0];
    return arg;
}

static vector<unsigned long> mkLens(vector<double> const &x,
				   vector<double> const &y,
				   vector<double> const &z)
{
    vector<unsigned long> arglen(3);
    arglen[0] = x.size();
    arglen[1] = y.size();
    arglen[2] = z.size();
    return arglen;
}

//Evaluate a VectorFunction that takes three arguments

vector<double>
veval(VectorFunction const *f, vector<double> const &x,
      vector<double> const &y, vector<double> const &z)
{
    return checkVEval(f, mkArgs(x,y,z), mkLens(x,y,z));
}

//Evaluate a VectorFunction that takes four arguments

static vector<double const *> mkArgs(vector<double> const &x,
				     vector<double> const &y,
				     vector<double> const &z,
				     vector<double> const &w)
{
    vector<double const *> arg(4);
    arg[0] = &x[0];
    arg[1] = &y[0];
    arg[2] = &z[0];
    arg[3] = &w[0];
    return arg;
}

static vector<unsigned long> mkLens(vector<double> const &x,
				    vector<double> const &y,
				    vector<double> const &z,
				    vector<double> const &w)
{
    vector<unsigned long> arglen(4);
    arglen[0] = x.size();
    arglen[1] = y.size();
    arglen[2] = z.size();
    arglen[3] = w.size();
    return arglen;
}

vector<double>
veval(VectorFunction const *f,
      vector<double> const &x, vector<double> const &y,
      vector<double> const &z, vector<double> const &w)
{
    return checkVEval(f, mkArgs(x,y,z,w), mkLens(x,y,z,w));
}

bool checkargs(VectorFunction const *f,
	       vector<double> const &x, vector<double> const &y,
	       vector<double> const &z, vector<double> const &w)
{
    return checkVArgs(f, mkArgs(x,y,z,w), mkLens(x,y,z,w));
}


/*
  Evaluate a VectorFunction that takes a single argument and returns a
  scalar
*/
double eval(VectorFunction const *f, vector<double> const &x)
{
    vector<double> ans = veval(f, x);
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.size());
    return ans[0];
}

/*
  Evaluate a VectorFunction that takes two arguments and returns a scalar
*/
double eval(VectorFunction const *f, vector<double> const &x, 
	    vector<double> const &y)
{
    vector<double> ans = veval(f, x, y);
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.size());
    return ans[0];
}

/*
  Evaluate a VectorFunction that takes three arguments and returns a scalar
*/
double eval(VectorFunction const *f, vector<double> const &x, 
	    vector<double> const &y, vector<double> const &z)
{
    vector<double> ans = veval(f, x, y, z);
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.size());
    return ans[0];
}


/*
  Array functions
*/

static vector<bool> discreteMask(vector<double const *> const &args,
				 vector<vector<unsigned long>> const &dims)
{
    vector<bool> out(args.size(), true);
    
    for (unsigned long i = 0; i < args.size(); ++i) {
	double const *v = args[i];
	unsigned long arglen = product(dims[i]);
	for (unsigned long j = 0; j < arglen; ++j) {
	    if (v[j] != floor(v[j] + 0.5)) {
		out[i] = false;
		break;
	    }
	}
    }

    return out;
}

static bool checkAArgs(ArrayFunction const *f,
		       vector<double const *> const &args,
		       vector<vector<unsigned long>> const &dims)
{
    return args.size() == dims.size() &&
	checkNPar(f, args.size())  &&
	f->checkParameterDim(dims) &&
	f->checkParameterDiscrete(discreteMask(args, dims)) &&
	f->checkParameterValue(args, dims);
}

pair<vector<double>, vector<unsigned long>>
checkAEval(ArrayFunction const *f,
	   vector<double const *> const &args,
	   vector<vector<unsigned long>> const &argdims)
{
    // Evaluate array function with checks
    CPPUNIT_ASSERT_MESSAGE(string("Valid arguments for ") + f->name(),
			   checkAArgs(f, args, argdims));
    vector<unsigned long> dim = f->dim(argdims, args);

    vector<double> value(product(dim));
    f->evaluate(value.data(), args, argdims);

    return pair<vector<double>, vector<unsigned long>> (value, dim);
}

static double const *getData(array_value const &a)
{
    return a.first.data();
}

static vector<unsigned long> const &getDims(array_value const &a)
{
    return a.second;
}

void addArgs(vector<double const *> const &v)
{
}

template<typename... Args>
void addArgs(vector<double const *> &v, array_value const &arg1, Args&... args)
{
    v.push_back(getData(arg1));
    addArgs(v, args...);
}

template<typename... Args>
vector<double const *> mkArgs(Args&... args)
{
    vector<double const *> v;
    addArgs(v, args...);
    return v;
}

void addDims(vector<vector<unsigned long>> const &v)
{
}

template<typename... Args>
void addDims(vector<vector<unsigned long>> &v, array_value const &arg1, Args&... args)
{
    v.push_back(getDims(arg1));
    addDims(v, args...);
}

template<typename... Args>
vector<vector<unsigned long>> mkDims(Args&... args)
{
    vector<vector<unsigned long>> v;
    addDims(v, args...);
    return v;
}

bool all_equal(array_value const &A, array_value const &B, double tol)
{
    if (A.second.size() != B.second.size()) return false;
    for (unsigned long j = 0; j < A.second.size(); ++j) {
	if (A.second[j] != B.second[j]) return false;
    }
    if (A.first.size() != B.first.size()) return false;
    for (unsigned long i = 0; i < A.first.size(); ++i) {
	if (abs(A.first[i] - B.first[i]) > tol) return false;
    }
    return true;
}

// Array function taking a single argument

array_value aeval(ArrayFunction const *f, array_value const &x)
{
    return checkAEval(f, mkArgs(x), mkDims(x));
}


bool checkargs(ArrayFunction const *f, array_value const &x)
{
    return checkAArgs(f, mkArgs(x), mkDims(x));
}

// Array function taking 2 arguments

array_value aeval(ArrayFunction const *f, array_value const &x,
		  array_value const &y)
{
    return checkAEval(f, mkArgs(x,y), mkDims(x,y));
}


bool checkargs(ArrayFunction const *f, array_value const &x,
	       array_value const &y)
{
    return checkAArgs(f, mkArgs(x,y), mkDims(x,y));
}

// Array function taking 3 arguments

array_value aeval(ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z)
{
    return checkAEval(f, mkArgs(x,y,z), mkDims(x,y,z));
}


bool checkargs(ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z)
{
    return checkAArgs(f, mkArgs(x,y,z), mkDims(x,y,z));
}

// Array function taking 4 arguments

array_value aeval(ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z,
		  array_value const &u)
{
    return checkAEval(f, mkArgs(x,y,z,u), mkDims(x,y,z,u));
}


bool checkargs(ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z,
	       array_value const &u)
{
    return checkAArgs(f, mkArgs(x,y,z,u), mkDims(x,y,z,u));
}

// Array function taking 5 arguments

array_value aeval(ArrayFunction const *f, array_value const &x,
		  array_value const &y, array_value const &z,
		  array_value const &u, array_value const &v)
{
    return checkAEval(f, mkArgs(x,y,z,u,v), mkDims(x,y,z,u,v));
}


bool checkargs(ArrayFunction const *f, array_value const &x,
	       array_value const &y, array_value const &z,
	       array_value const &u, array_value const &v)
{
    return checkAArgs(f, mkArgs(x,y,z,u,v), mkDims(x,y,z,u,v));
}

/* Array functions returning a scalar */

double eval(ArrayFunction const *f, array_value const &x)
{
    array_value ans = checkAEval(f, mkArgs(x), mkDims(x));
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);
    return ans.first[0];
}

double eval(ArrayFunction const *f, array_value const &x,
	    array_value const &y)
{
    array_value ans = checkAEval(f, mkArgs(x,y), mkDims(x,y));
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);
    return ans.first[0];
}

double eval(ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z)
{
    array_value ans = checkAEval(f, mkArgs(x,y,z), mkDims(x,y,z));
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);
    return ans.first[0];
}

double eval(ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z,
	    array_value const &u)
{
    array_value ans = checkAEval(f, mkArgs(x,y,z,u), mkDims(x,y,z,u));
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);
    return ans.first[0];
}

double eval(ArrayFunction const *f, array_value const &x,
	    array_value const &y, array_value const &z,
	    array_value const &u, array_value const &v)
{
    array_value ans = checkAEval(f, mkArgs(x,y,z,u,v), mkDims(x,y,z,u,v));
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);
    return ans.first[0];
}
