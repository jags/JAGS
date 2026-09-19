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

bool checkArgs(ScalarFunction const *f,
	       vector<double const *> const &args)
{
    //Check that arguments are valid
    return checkNPar(f, args.size()) &&
	f->checkParameterDiscrete(discreteMask(args)) &&
	f->checkParameterValue(args);
}

double Eval(ScalarFunction const *f,
	    vector<double const *> const &args)
{
    //Evaluate scalar function with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args));
    return f->evaluate(args);
}

double Gradient(ScalarFunction const *f,
		vector<double const *> const &args,
		unsigned long i)
{
    //Evaluate gradient with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));
    return f->gradient(args, i);
}

double NumGradient(ScalarFunction const *f,
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
    double y1 = Eval(f, args1);
    args0[i] = *args[i] + delta;
    double y2 = Eval(f, args1);
	
    return (y2 - y1)/(2*delta);
}

/* Vector functions */

bool all_equal(vector<double> const &u, vector<double> const &v, double tol)
{
    if (u.size() != v.size()) return false;
    for (unsigned long i = 0; i < u.size(); ++i) {
	if (abs(u[i] - v[i]) > tol) return false;
    }
    return true;
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
			   
bool checkArgs(VectorFunction const *f,
		vector<double const *> const &args,
		vector<unsigned long> const &arglen)
{
    return args.size() == arglen.size() &&
	checkNPar(f, args.size()) &&
	f->checkParameterLength(arglen) &&
	f->checkParameterDiscrete(discreteMask(args, arglen)) &&
	f->checkParameterValue(args, arglen);
}

vector<double> VEval(VectorFunction const *f,
		     vector<double const *> const &args,
		     vector<unsigned long> const &arglen)
{
    // Evaluate vector function with checks
    CPPUNIT_ASSERT_MESSAGE(string("Valid arguments for ") + f->name(),
			   checkArgs(f, args, arglen));
    vector<double> ans(f->length(arglen, args));
    f->evaluate(&ans[0], args, arglen);
    return ans;
}

double Eval(VectorFunction const *f,
	    vector<double const *> const &args,
	    vector<unsigned long> const &arglen)
{
    vector<double> ans = VEval(f, args, arglen);
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.size());
    return ans[0];
}

vector<double> VGradient(VectorFunction const *f,
			 vector<double const *> const &args,
			 vector<unsigned long> const &arglen,
			 unsigned long i)
{
    //Evaluate gradient with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args, arglen));
    unsigned long n = f->length(arglen, args);
    unsigned long m = arglen[i];
    vector<double> ans(n * m, 0);

    f->gradient(ans.data(), args, arglen, i);
    return ans;
}

vector<double> VNumGradient(VectorFunction const *f,
			    vector<double const *> const &args,
			    vector<unsigned long> const &arglen,
			    unsigned long i, double delta)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args, arglen));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));

    //Dimensions of answer matrix
    unsigned long n = f->length(arglen, args);
    unsigned long m = arglen[i];

    //Create mutable copy of the ith argument
    vector<double const *> args1(args);
    vector<double> argi(m);
    copy(args[i], args[i] + m, argi.begin());
    args1[i] = argi.data();


    vector<double> ans(n * m, 0);
    for (unsigned long j = 0; j < m; ++j) {
	argi[j] = args[i][j] - delta;
	vector<double> y1 = VEval(f, args1, arglen);
	argi[j] = args[i][j] + delta;
	vector<double> y2 = VEval(f, args1, arglen);
	argi[j] = args[i][j];
	for (unsigned long k = 0; k < n; ++k) {
	    ans[j*n + k] = (y2[k] - y1[k])/(2*delta);
	}
    }

    return ans;
}

/*
  Array functions
*/

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

bool checkArgs(ArrayFunction const *f,
	       vector<double const *> const &args,
	       vector<vector<unsigned long>> const &dims)
{
    return args.size() == dims.size() &&
	checkNPar(f, args.size())  &&
	f->checkParameterDim(dims) &&
	f->checkParameterDiscrete(discreteMask(args, dims)) &&
	f->checkParameterValue(args, dims);
}

array_value AEval(ArrayFunction const *f,
		  vector<double const *> const &args,
		  vector<vector<unsigned long>> const &argdims)
{
    // Evaluate array function with checks
    CPPUNIT_ASSERT_MESSAGE(string("Valid arguments for ") + f->name(),
			   checkArgs(f, args, argdims));
    vector<unsigned long> dim = f->dim(argdims, args);

    vector<double> value(product(dim));
    f->evaluate(value.data(), args, argdims);

    return array_value(value, dim);
}

double Eval(ArrayFunction const *f,
	    vector<double const *> const &args,
	    vector<vector<unsigned long>> const &argdims)
{
    array_value ans = AEval(f, args, argdims);
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.first.size());
    CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), ans.second.size());
    CPPUNIT_ASSERT_EQUAL(1UL, ans.second[0]);

    return ans.first[0];
}

vector<double> VGradient(ArrayFunction const *f,
			 vector<double const *> const &args,
			 vector<vector<unsigned long>> const &argdims,
			 unsigned long i)
{
    //Evaluate gradient with checks
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args, argdims));
    unsigned long n = product(f->dim(argdims, args));
    unsigned long m = product(argdims[i]);
    vector<double> ans(n * m, 0);

    f->gradient(ans.data(), args, argdims, i);
    return ans;
}

vector<double> VNumGradient(ArrayFunction const *f,
			    vector<double const *> const &args,
			    vector<vector<unsigned long>> const &argdims,
			    unsigned long i, double delta)
{
    CPPUNIT_ASSERT_MESSAGE(f->name(), checkArgs(f, args, argdims));
    CPPUNIT_ASSERT_MESSAGE(f->name(), f->hasGradient(i));

    //Dimensions of answer matrix
    unsigned long n = product(f->dim(argdims, args));
    unsigned long m = product(argdims[i]);

    //Create mutable copy of the ith argument
    vector<double const *> args1(args);
    vector<double> argi(m);
    copy(args[i], args[i] + m, argi.begin());
    args1[i] = argi.data();


    vector<double> ans(n * m, 0);
    for (unsigned long j = 0; j < m; ++j) {
	argi[j] = args[i][j] - delta;
	vector<double> y1 = AEval(f, args1, argdims).first;
	argi[j] = args[i][j] + delta;
	vector<double> y2 = AEval(f, args1, argdims).first;
	argi[j] = args[i][j];
	for (unsigned long k = 0; k < n; ++k) {
	    ans[j*n + k] = (y2[k] - y1[k])/(2*delta);
	}
    }
    
    return ans;
}
