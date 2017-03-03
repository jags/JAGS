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
    CPPUNIT_TEST_SUITE_END(  );

    jags::RNG *_rng;

    jags::RScalarDist *_dbetabin;

    void rscalar_rpq(jags::RScalarDist const *dist, 
		     std::vector<double const *> const &par);

    void dkwtest(jags::RScalarDist const *dist,
		 std::vector<double const *> const &par,
		 unsigned int N=10000, double pthresh=0.001);
    
  public:
    void setUp();
    void tearDown();

    void npar();
    void name();
    void alias();
    void rscalar();

    void dkw();
};

#endif /* MIX_DIST_TEST_H */
