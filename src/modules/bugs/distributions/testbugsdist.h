#ifndef BUGS_DIST_TEST_H
#define BUGS_DIST_TEST_H

namespace jags {
    class RScalarDist;
    class ScalarDist;
    class VectorDist;
    class ArrayDist;
    struct RNG;
}

#include <cppunit/extensions/HelperMacros.h>
#include <testlib.h>

class BugsDistTest : public CppUnit::TestFixture, public JAGSFixture
{

    CPPUNIT_TEST_SUITE( BugsDistTest );
    CPPUNIT_TEST( npar );
    CPPUNIT_TEST( name );
    CPPUNIT_TEST( alias );
    CPPUNIT_TEST( rscalar );
    CPPUNIT_TEST( kl );
    CPPUNIT_TEST( dkw );
    CPPUNIT_TEST_SUITE_END(  );

    jags::RNG *_rng;

    jags::ScalarDist *_dbern;
    jags::RScalarDist *_dbeta;
    jags::RScalarDist *_dbin;
    jags::VectorDist *_dcat;
    jags::RScalarDist *_dchisqr;
    jags::VectorDist *_ddirch;
    jags::RScalarDist *_ddexp;
    jags::RScalarDist *_dexp;
    jags::RScalarDist *_df;
    jags::RScalarDist *_dgamma;
    jags::RScalarDist *_dgengamma;
    jags::RScalarDist *_dhyper;
    jags::VectorDist *_dinterval;    
    jags::RScalarDist *_dlnorm;
    jags::RScalarDist *_dlogis;
    jags::ArrayDist *_dmnorm;
    jags::ArrayDist *_dmnormvc;
    jags::ArrayDist *_dmt;
    jags::VectorDist *_dmulti;
    jags::RScalarDist *_dnchisqr;
    jags::RScalarDist *_dnegbin;
    jags::RScalarDist *_dgampois;
    jags::RScalarDist *_dnorm;
    jags::RScalarDist *_dnt;
    jags::RScalarDist *_dpar;
    jags::RScalarDist *_dpois;
    jags::ScalarDist *_dround;
    jags::VectorDist *_dsample;
    jags::ArrayDist *_dsum;
    jags::RScalarDist *_dt;
    jags::ScalarDist *_dunif;
    jags::RScalarDist *_dweib;
    jags::ArrayDist *_dwish;

    void rscalar_rpq(jags::RScalarDist const *dist, 
		     std::vector<double const *> const &par);

    void rscalar_trunclik(jags::RScalarDist const *dist, 
			  std::vector<double const *> const &par);

    void kl_scalar(jags::ScalarDist const *dist, 
		   std::vector<double const *> const &par0,
		   std::vector<double const *> const &par1);

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

    void kl();
    void dkw();
};

#endif /* BUGS_DIST_TEST_H */
