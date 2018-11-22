#include <function/testfun.h>
#include <util/nainf.h>
#include <util/logical.h>

#include <cppunit/extensions/HelperMacros.h>

using jags::ScalarFunction;
using jags::VectorFunction;
using jags::Function;
using jags::anyTrue;
using jags::allTrue;

#include <climits>
#include <cmath>
#include <algorithm>

using std::vector;
using std::string;
using std::copy;
using std::floor;

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

    CPPUNIT_ASSERT_MESSAGE(f->name(), !jags_isnan(eval(f, lower)));
    CPPUNIT_ASSERT_MESSAGE(f->name(), !jags_isnan(eval(f, upper)));

    if (jags_finite(upper)) {
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
    if (jags_finite(lower)) {
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


/* Evaluate a VectorFunction that takes a single argument */
static vector<double const *> mkArgs(vector<double> const &x)
{
    return vector<double const *>(1, &x[0]);
    vector<unsigned long> arglen(1, x.size());
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

