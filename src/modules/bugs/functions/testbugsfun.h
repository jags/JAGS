#ifndef BUGS_FUN_TEST_H
#define BUGS_FUN_TEST_H

namespace jags {
    class ScalarFunction;
    class VectorFunction;
    class ArrayFunction;
    class LinkFunction;
}

#include <cppunit/extensions/HelperMacros.h>
#include <testlib.h>
#include <function/testfun.h>

class BugsFunTest : public CppUnit::TestFixture, public JAGSFixture
{
    //FIXME: Failing tests are commented out
    CPPUNIT_TEST_SUITE( BugsFunTest );
    CPPUNIT_TEST( npar );
    CPPUNIT_TEST( name );
    CPPUNIT_TEST( alias );
    CPPUNIT_TEST( trig );

    CPPUNIT_TEST( hyper );
    CPPUNIT_TEST( link );
    //CPPUNIT_TEST( summary );
    //CPPUNIT_TEST( math );
    CPPUNIT_TEST( lossy );
    CPPUNIT_TEST( slap );

    CPPUNIT_TEST( additive );
    CPPUNIT_TEST( linear );
    CPPUNIT_TEST( scale );
    //CPPUNIT_TEST( sort );
    CPPUNIT_TEST( matrix );

    //CPPUNIT_TEST( inprod );
    CPPUNIT_TEST( ifelse );
    CPPUNIT_TEST( interplin );
    CPPUNIT_TEST( discrete );
    CPPUNIT_TEST( combine );
    CPPUNIT_TEST( rep );
    CPPUNIT_TEST_SUITE_END();

    jags::ScalarFunction *_abs;
    jags::ScalarFunction *_arccos;    
    jags::ScalarFunction *_arccosh;
    jags::ScalarFunction *_arcsin;
    jags::ScalarFunction *_arcsinh;
    jags::ScalarFunction *_arctan;
    jags::ScalarFunction *_arctanh;
    jags::ScalarFunction *_cloglog;
    jags::VectorFunction *_combine;
    jags::ScalarFunction *_cos;
    jags::ScalarFunction *_cosh;
    jags::VectorFunction *_dinterval;
    jags::ScalarFunction *_dround;
    jags::ArrayFunction  *_dsum;
    jags::LinkFunction   *_exp;
    jags::LinkFunction   *_icloglog;
    jags::ScalarFunction *_ifelse;
    jags::LinkFunction   *_ilogit;
    jags::VectorFunction *_inprod;
    jags::VectorFunction *_interplin;
    jags::ArrayFunction  *_inverse;
    jags::ArrayFunction  *_logdet;
    jags::ScalarFunction *_logfact;
    jags::ScalarFunction *_loggam;
    jags::ScalarFunction *_log;
    jags::ScalarFunction *_logit;
    jags::ArrayFunction  *_matmult;
    jags::VectorFunction *_max;
    jags::VectorFunction *_mean;
    jags::VectorFunction *_min;
    jags::VectorFunction *_order;
    jags::LinkFunction   *_phi;
    jags::ScalarFunction *_probit;
    jags::VectorFunction *_prod;
    jags::VectorFunction *_rank;
    jags::VectorFunction *_rep;
    jags::ScalarFunction *_round;
    jags::VectorFunction *_sd;
    jags::ScalarFunction *_sin;
    jags::ScalarFunction *_sinh;
    jags::VectorFunction *_sort;
    jags::ScalarFunction *_sqrt;
    jags::ScalarFunction *_step;
    jags::VectorFunction *_sum;
    jags::ScalarFunction *_tan;
    jags::ScalarFunction *_tanh;
    jags::ArrayFunction  *_transpose;
    jags::ScalarFunction *_trunc;
    jags::ScalarFunction *_yytrunc;

    void link(jags::ScalarFunction const *f, jags::LinkFunction const *link,
	      double lower, double upper, int N);


    /*
    void summary(double const *v, unsigned int N);
    void summary(double const *v1, unsigned int N1,
		 double const *v2, unsigned int N2);
    */

    void summary(std::vector<double> const &v);
    void summary(std::vector<double> const &v1, std::vector<double> const &v2);


    /*
    template<size_t N>
    void summary(double const (&x)[N]) {
	summary(x, N);
    }

    template<size_t N1, size_t N2>
    void summary(double const (&x)[N1], double const (&y)[N2]) {
	summary(x, N1, y, N2);
    }
    */
    
    template<size_t N>
    void summary(double const (&x)[N]) {
	summary(mkVec(x));
    }

    template<size_t N1, size_t N2>
    void summary(double const (&x)[N1], double const (&y)[N2]) {
	summary(mkVec(x), mkVec(y));
    }
    
    void trig(double v);
    void hyper(double v);

  public:
    void setUp();
    void tearDown();

    void npar();
    void name();
    void alias();
    void link();
    void summary();
    void trig();
    void hyper();
    void math();
    void lossy();
    void slap();
    void additive();
    void linear();
    void scale();
    void sort();
    void power();
    void matrix();
    void inprod();
    void ifelse();
    void discrete();
    void interplin();
    void combine();
    void rep();
};

#endif  // BUGS_FUN_TEST_H


