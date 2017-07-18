#include "testbasefun.h"

#include "Add.h"
#include "And.h"
#include "Divide.h"
#include "Equal.h"
#include "GreaterOrEqual.h"
#include "GreaterThan.h"
#include "LessOrEqual.h"
#include "LessThan.h"
#include "Multiply.h"
#include "Neg.h"
#include "NotEqual.h"
#include "Not.h"
#include "Or.h"
#include "Pow.h"
#include "Seq.h"
#include "Subtract.h"

#include <function/testfun.h>
#include <util/nainf.h>
#include <util/integer.h>
using jags::checkInteger;

#include <cmath>
using std::pow;
using std::vector;
using std::string;

#include <climits>

void BaseFunTest::setUp()
{
    _add = new jags::base::Add;
    _and = new jags::base::And;
    _divide = new jags::base::Divide;
    _equal = new jags::base::Equal;
    _geq = new jags::base::GreaterOrEqual;
    _gt = new jags::base::GreaterThan;
    _leq = new jags::base::LessOrEqual;
    _lt = new jags::base::LessThan;
    _multiply = new jags::base::Multiply;
    _neg = new jags::base::Neg;
    _neq = new jags::base::NotEqual;
    _not = new jags::base::Not;
    _or = new jags::base::Or;
    _pow = new jags::base::Pow;
    _seq = new jags::base::Seq;
    _subtract = new jags::base::Subtract;
}

void BaseFunTest::tearDown()
{
    delete _add;
    delete _and;
    delete _divide;
    delete _equal;
    delete _geq;
    delete _gt;
    delete _leq;
    delete _lt;
    delete _multiply;
    delete _neg;
    delete _neq;
    delete _not;
    delete _or;
    delete _pow;
    delete _seq;
    delete _subtract;
}

void BaseFunTest::name()
{
    CPPUNIT_ASSERT_EQUAL(string("+"), _add->name());
    CPPUNIT_ASSERT_EQUAL(string("&&"), _and->name());
    CPPUNIT_ASSERT_EQUAL(string("/"), _divide->name());
    CPPUNIT_ASSERT_EQUAL(string("=="), _equal->name());
    CPPUNIT_ASSERT_EQUAL(string(">="), _geq->name());
    CPPUNIT_ASSERT_EQUAL(string(">"), _gt->name());
    CPPUNIT_ASSERT_EQUAL(string("<="), _leq->name());
    CPPUNIT_ASSERT_EQUAL(string("<"), _lt->name());
    CPPUNIT_ASSERT_EQUAL(string("*"), _multiply->name());
    //Unary minus is converted to the token NEG to avoid
    //confusing with binary minus
    CPPUNIT_ASSERT_EQUAL(string("NEG"), _neg->name());
    CPPUNIT_ASSERT_EQUAL(string("!="), _neq->name());
    CPPUNIT_ASSERT_EQUAL(string("!"), _not->name());
    CPPUNIT_ASSERT_EQUAL(string("||"), _or->name());
    CPPUNIT_ASSERT_EQUAL(string("^"), _pow->name());
    CPPUNIT_ASSERT_EQUAL(string("-"), _subtract->name());
    CPPUNIT_ASSERT_EQUAL(string(":"), _seq->name());
}

void BaseFunTest::alias()
{
    //None of the built-in functions has an alias, except pow
    CPPUNIT_ASSERT_EQUAL(string(""), _add->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _and->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _divide->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _equal->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _geq->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _gt->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _leq->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _lt->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _multiply->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _neg->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _neq->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _not->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _or->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _subtract->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _seq->alias());
    CPPUNIT_ASSERT_EQUAL(string("pow"), _pow->alias());
}

void BaseFunTest::arithmetic1(double v)
{
    //Arithmetic functions taking a single argument
    CPPUNIT_ASSERT_EQUAL(_neg->npar(), 1U);
    CPPUNIT_ASSERT(checkNPar(_add, 1));
    CPPUNIT_ASSERT(checkNPar(_multiply, 1));

    CPPUNIT_ASSERT_EQUAL(-v, eval(_neg, v));
    CPPUNIT_ASSERT_EQUAL(v, eval(_add, v));
    CPPUNIT_ASSERT_EQUAL(v, eval(_multiply, v));
}

void BaseFunTest::arithmetic2(double v1, double v2)
{
    //Arithmetic functions taking two arguments
    CPPUNIT_ASSERT(checkNPar(_add, 2));
    CPPUNIT_ASSERT(checkNPar(_multiply, 2));
    CPPUNIT_ASSERT_EQUAL(_subtract->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_divide->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_pow->npar(), 2U);

    CPPUNIT_ASSERT_EQUAL(v1 + v2, eval(_add, v1, v2));
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v1 - v2, eval(_subtract,  v1, v2), tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v1 * v2, eval(_multiply, v1, v2), tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v1 / v2, eval(_divide, v1, v2), tol);
    if (v1 < 0 && !checkInteger(v2)) {
	CPPUNIT_ASSERT_ASSERTION_FAIL(eval(_pow, v1, v2));
    }
    else {
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pow(v1, v2), eval(_pow, v1, v2), tol);
    }
}

void BaseFunTest::arithmetic3(double v1, double v2, double v3)
{
    //Arithmetic functions taking 3 arguments
    CPPUNIT_ASSERT(checkNPar(_add, 3));
    CPPUNIT_ASSERT(checkNPar(_multiply, 3));

    vector<double const *> args(3);
    args[0] = &v1;
    args[1] = &v2;
    args[2] = &v3;

    CPPUNIT_ASSERT_EQUAL(_add->evaluate(args), v1 + v2 + v3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_multiply->evaluate(args), v1 * v2 * v3, tol);
}

void BaseFunTest::arithmetic()
{
    arithmetic1(0.0);
    arithmetic1(-3.7);
    arithmetic1(28.16);

    arithmetic2( 1.0,  2.0);
    arithmetic2(-1.0,  3.0);
    arithmetic2(-2.5,  8.3);
    arithmetic2(-3.7, -1.9);

    arithmetic3( 0.0, 0.0,    0.0);
    arithmetic3( 1.0, 2.0,    3.0);
    arithmetic3(-1.0, 3.0,  -19.2);
    arithmetic3(-2.5,  8.3, -8.62);
    arithmetic3(-3.7, -1.9, -3.8);
}

void BaseFunTest::logical1()
{
    //Logical Not takes one argument
    CPPUNIT_ASSERT_EQUAL(_not->npar(), 1U);

    CPPUNIT_ASSERT_EQUAL(1.0, eval(_not, 0));
    CPPUNIT_ASSERT_EQUAL(0.0, eval(_not, 1));
    //Any non-zero value should be interpreted as TRUE
    CPPUNIT_ASSERT_EQUAL(0.0, eval(_not, 2));
    CPPUNIT_ASSERT_EQUAL(0.0, eval(_not, -0.5));
}

void BaseFunTest::logical2()
{
    //Logical And and Logical Or take two arguments
    CPPUNIT_ASSERT_EQUAL(_and->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_or->npar(), 2U);

    double TRUE=1, FALSE=0;

    CPPUNIT_ASSERT_EQUAL(FALSE, eval(_and, FALSE, FALSE));
    CPPUNIT_ASSERT_EQUAL(FALSE, eval(_and, FALSE, TRUE));
    CPPUNIT_ASSERT_EQUAL(FALSE, eval(_and, TRUE, FALSE));
    CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_and, TRUE, TRUE));

    CPPUNIT_ASSERT_EQUAL(FALSE, eval(_or, FALSE, FALSE));
    CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, FALSE, TRUE));
    CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, TRUE, FALSE));
    CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, TRUE, TRUE));

    //Any non-zero value should be interpreted as "true".
    //Try it with a few arbitrary non-zero values.
    double truevals[4] = {-3.872, 1e100, DBL_EPSILON, 3};

    for (unsigned int i = 0; i < 4; ++i) {

	CPPUNIT_ASSERT_EQUAL(FALSE, eval(_and, FALSE, truevals[i]));
	CPPUNIT_ASSERT_EQUAL(FALSE, eval(_and, truevals[i], FALSE));
	CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_and, truevals[i], truevals[i]));
	
	CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, FALSE, truevals[i]));
	CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, truevals[i], FALSE));
	CPPUNIT_ASSERT_EQUAL(TRUE,  eval(_or, truevals[i], truevals[i]));
    }
}

void BaseFunTest::logical()
{
    logical1();
    logical2();
}

void BaseFunTest::comparison2(double v1, double v2)
{
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 == v2), eval(_equal, v1, v2));
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 >= v2), eval(_geq, v1, v2));
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 > v2),  eval(_gt, v1, v2));
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 <= v2), eval(_leq, v1, v2));
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 < v2),  eval(_lt, v1, v2));
    CPPUNIT_ASSERT_EQUAL(static_cast<double>(v1 != v2), eval(_neq, v1, v2));
}

void BaseFunTest::comparison()
{
    //All of the comparison functions take 2 arguments
    CPPUNIT_ASSERT_EQUAL(_equal->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_geq->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_gt->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_leq->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_lt->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_neq->npar(), 2U);

    comparison2( 0.0,  0.0);
    comparison2( 1.0,  2.0);
    comparison2( 1.0,  1.0);
    comparison2(-1.0,  3.0);
    comparison2(-2.5,  8.3);
    comparison2(-3.7, -1.9);
    comparison2( -2.8, -2.8);
    comparison2( 0.0, DBL_EPSILON);
    comparison2( DBL_EPSILON, 0.0);
}

void BaseFunTest::discrete()
{

    CPPUNIT_ASSERT(isdiscrete(_not, 1, always));
    CPPUNIT_ASSERT(isdiscrete(_and, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_equal, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_geq, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_gt, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_leq, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_lt, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_neq, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_not, 1, always));
    CPPUNIT_ASSERT(isdiscrete(_or, 2, always));

    CPPUNIT_ASSERT(isdiscrete(_multiply, 3, all));
    CPPUNIT_ASSERT(isdiscrete(_add, 3, all));
    CPPUNIT_ASSERT(isdiscrete(_divide, 2, never));
    CPPUNIT_ASSERT(isdiscrete(_subtract, 2, all));

    CPPUNIT_ASSERT(isdiscrete(_seq, 2, always));
    CPPUNIT_ASSERT(isdiscrete(_pow, 2, never));

}

void BaseFunTest::slp()
{
    //Functions that are not scale, linear or power functions
    CPPUNIT_ASSERT(neverclosed(_and, 2));
    CPPUNIT_ASSERT(neverclosed(_equal, 2));
    CPPUNIT_ASSERT(neverclosed(_geq, 2));
    CPPUNIT_ASSERT(neverclosed(_gt, 2));
    CPPUNIT_ASSERT(neverclosed(_leq, 2));
    CPPUNIT_ASSERT(neverclosed(_gt, 2));
    CPPUNIT_ASSERT(neverclosed(_leq, 2));
    CPPUNIT_ASSERT(neverclosed(_lt, 2));
    CPPUNIT_ASSERT(neverclosed(_neq, 2));
    CPPUNIT_ASSERT(neverclosed(_not, 1));
    CPPUNIT_ASSERT(neverclosed(_or, 2));
    CPPUNIT_ASSERT(neverclosed(_seq, 2));
}

void BaseFunTest::linear()
{
    CPPUNIT_ASSERT(_add->isLinear(TT, vector<bool>()));
    CPPUNIT_ASSERT(_add->isLinear(TF, vector<bool>()));
    CPPUNIT_ASSERT(_add->isLinear(FT, vector<bool>()));
    CPPUNIT_ASSERT(_add->isLinear(TF, FT));
    CPPUNIT_ASSERT(_add->isLinear(TF, FF));
    CPPUNIT_ASSERT(_add->isLinear(FT, TF));
    CPPUNIT_ASSERT(_add->isLinear(FT, FF));

    CPPUNIT_ASSERT(_subtract->isLinear(TT, vector<bool>()));
    CPPUNIT_ASSERT(_subtract->isLinear(TF, vector<bool>()));
    CPPUNIT_ASSERT(_subtract->isLinear(FT, vector<bool>()));
    CPPUNIT_ASSERT(_subtract->isLinear(TF, FT));
    CPPUNIT_ASSERT(_subtract->isLinear(TF, FF));
    CPPUNIT_ASSERT(_subtract->isLinear(FT, TF));
    CPPUNIT_ASSERT(_subtract->isLinear(FT, FF));

    CPPUNIT_ASSERT(!_divide->isLinear(TT, vector<bool>()));
    CPPUNIT_ASSERT(!_divide->isLinear(FT, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isLinear(TF, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isLinear(TF, FT));

    CPPUNIT_ASSERT(!_multiply->isLinear(TT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isLinear(FT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isLinear(TF, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isLinear(TF, FT));
    CPPUNIT_ASSERT(_multiply->isLinear(FT, TF));
    
    CPPUNIT_ASSERT(neverlinear(_pow, 2));
}

void BaseFunTest::scale()
{
    CPPUNIT_ASSERT(_add->isScale(TT, vector<bool>()));
    CPPUNIT_ASSERT(!_add->isScale(TF, vector<bool>()));
    CPPUNIT_ASSERT(!_add->isScale(FT, vector<bool>()));

    CPPUNIT_ASSERT(_subtract->isScale(TT, vector<bool>()));
    CPPUNIT_ASSERT(!_subtract->isScale(TF, vector<bool>()));
    CPPUNIT_ASSERT(!_subtract->isScale(FT, vector<bool>()));

    CPPUNIT_ASSERT(!_divide->isScale(TT, vector<bool>()));
    CPPUNIT_ASSERT(!_divide->isScale(FT, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isScale(TF, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isScale(TF, FT));

    CPPUNIT_ASSERT(!_multiply->isScale(TT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isScale(FT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isScale(TF, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isScale(TF, FT));
    CPPUNIT_ASSERT(_multiply->isScale(FT, TF));

    CPPUNIT_ASSERT(neverscale(_pow, 2));
}

void BaseFunTest::power()
{
    CPPUNIT_ASSERT(neverpow(_add, 2));
    CPPUNIT_ASSERT(neverpow(_add, 3));
    CPPUNIT_ASSERT(neverpow(_subtract, 2));

    CPPUNIT_ASSERT(_divide->isPower(TT, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isPower(FT, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isPower(FT, TF));
    CPPUNIT_ASSERT(_divide->isPower(TF, vector<bool>()));
    CPPUNIT_ASSERT(_divide->isPower(TF, FT));

    CPPUNIT_ASSERT(_multiply->isPower(TT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isPower(FT, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isPower(FT, TF));
    CPPUNIT_ASSERT(_multiply->isPower(TF, vector<bool>()));
    CPPUNIT_ASSERT(_multiply->isPower(TF, FT));

    CPPUNIT_ASSERT(_pow->isPower(TF, vector<bool>()));
    CPPUNIT_ASSERT(_pow->isPower(TF, FT));
    CPPUNIT_ASSERT(!_pow->isPower(FT, vector<bool>()));
}

void BaseFunTest::seq()
{
    for (int i = -5; i < 6; ++i) {
	for (int j = -6; j < 5; ++j) {
	    int len = (j < i) ? 0 : j - i + 1;
	    vector<double> out = veval(_seq, i, j);
	    CPPUNIT_ASSERT(out.size() == static_cast<unsigned int>(len));
	    for (int k = 0; k < len; ++k) {
		CPPUNIT_ASSERT_EQUAL(out[k], static_cast<double>(i + k));
	    }
	}
    }
				     
}
