#ifndef MIX_DIST_TEST_H
#define MIX_DIST_TEST_H

namespace jags {
    class RScalarDist;
    class ScalarDist;
    class VectorDist;
    class ArrayDist;
    struct RNG;
}

#include <cppunit/extensions/HelperMacros.h>
#include <testlib.h>

class MixDistTest : public CppUnit::TestFixture, public JAGSFixture
{

    CPPUNIT_TEST_SUITE( MixDistTest );
    CPPUNIT_TEST( npar );
    CPPUNIT_TEST( name );
    CPPUNIT_TEST( alias );
    CPPUNIT_TEST( rscalar );
    CPPUNIT_TEST( dkw );
    CPPUNIT_TEST( normmix );
    CPPUNIT_TEST_SUITE_END(  );

    jags::RNG *_rng;

    jags::RScalarDist *_dbetabin;
    jags::VectorDist *_dnormmix;

    void rscalar_rpq(jags::RScalarDist const *dist, 
		     std::vector<double const *> const &par);

    void dkwtest(jags::RScalarDist const *dist,
		 std::vector<double const *> const &par,
		 unsigned int N=10000, double pthresh=0.001);

    void test_mean_normmix(std::vector<double const *> const &par,
			   std::vector<unsigned long> const &lengths,
			   unsigned int N);

    void test_var_normmix(std::vector<double const *> const &par,
			  std::vector<unsigned long> const &len,
			  unsigned int N);

  public:
    void setUp();
    void tearDown();

    void npar();
    void name();
    void alias();
    void rscalar();
    
    void dkw();

    void normmix();
};

#endif /* MIX_DIST_TEST_H */
