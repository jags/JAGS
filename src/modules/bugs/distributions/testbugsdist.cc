#include "testbugsdist.h"

#include "DBern.h"
#include "DBeta.h"
#include "DBin.h"
#include "DCat.h"
#include "DChisqr.h"
#include "DDexp.h"
#include "DDirch.h"
#include "DExp.h"
#include "DF.h"
#include "DGamma.h"
#include "DGenGamma.h"
#include "DHyper.h"
#include "DInterval.h"
#include "DLnorm.h"
#include "DLogis.h"
#include "DMNorm.h"
#include "DMNormVC.h"
#include "DMT.h"
#include "DMulti.h"
#include "DNChisqr.h"
#include "DNegBin.h"
#include "DGamPois.h"
#include "DNorm.h"
#include "DNT.h"
#include "DPar.h"
#include "DPois.h"
#include "DRound.h"
#include "DSample.h"
#include "DSum.h"
#include "DT.h"
#include "DUnif.h"
#include "DWeib.h"
#include "DWish.h"

#include <MersenneTwisterRNG.h>
#include <util/nainf.h>
#include <JRmath.h>

#include <cmath>
#include <set>
#include <sstream>
#include <algorithm>

using std::string;
using std::vector;
using std::sqrt;
using std::max;
using std::min;
using std::multiset;
using std::abs;
using std::ostringstream;
using std::sort;

using jags::ScalarDist;
using jags::RScalarDist;

void BugsDistTest::setUp() {

    _rng = new jags::base::MersenneTwisterRNG(1234567, 
					      jags::KINDERMAN_RAMAGE);

    _dbern = new jags::bugs::DBern();
    _dbeta = new jags::bugs::DBeta();
    _dbin = new jags::bugs::DBin();
    _dcat = new jags::bugs::DCat();
    _dchisqr = new jags::bugs::DChisqr();
    _ddirch = new jags::bugs::DDirch();
    _ddexp = new jags::bugs::DDexp();
    _dexp = new jags::bugs::DExp();
    _df = new jags::bugs::DF();
    _dgamma = new jags::bugs::DGamma();
    _dgengamma = new jags::bugs::DGenGamma();
    _dhyper = new jags::bugs::DHyper();
    _dinterval = new jags::bugs::DInterval();    
    _dlnorm = new jags::bugs::DLnorm();
    _dlogis = new jags::bugs::DLogis();
    _dmnorm = new jags::bugs::DMNorm();
    _dmnormvc = new jags::bugs::DMNormVC();
    _dmt = new jags::bugs::DMT();
    _dmulti = new jags::bugs::DMulti();
    _dnchisqr = new jags::bugs::DNChisqr();
    _dnegbin = new jags::bugs::DNegBin();
    _dgampois = new jags::bugs::DGamPois();
    _dnorm = new jags::bugs::DNorm();
    _dnt = new jags::bugs::DNT();
    _dpar = new jags::bugs::DPar();
    _dpois = new jags::bugs::DPois();
    _dround = new jags::bugs::DRound();
    _dsample = new jags::bugs::DSample();
    _dsum = new jags::bugs::DSum();
    _dt = new jags::bugs::DT();
    _dunif = new jags::bugs::DUnif();
    _dweib = new jags::bugs::DWeib();
    _dwish = new jags::bugs::DWish();

}

void BugsDistTest::tearDown() {

    delete _rng;

    delete _dbern;
    delete _dbeta;
    delete _dbin;
    delete _dcat;
    delete _dchisqr;
    delete _ddirch;
    delete _ddexp;
    delete _dexp;
    delete _df;
    delete _dgamma;
    delete _dgengamma;
    delete _dhyper;
    delete _dinterval;    
    delete _dlnorm;
    delete _dlogis;
    delete _dmnorm;
    delete _dmnormvc;
    delete _dmt;
    delete _dmulti;
    delete _dnchisqr;
    delete _dnegbin;
    delete _dgampois;
    delete _dnorm;
    delete _dnt;
    delete _dpar;
    delete _dpois;
    delete _dround;
    delete _dsample;
    delete _dsum;
    delete _dt;
    delete _dunif;
    delete _dweib;
    delete _dwish;
}

void BugsDistTest::npar()
{
    //CPPUNIT_ASSERT_MESSAGE("npar check", false);

    CPPUNIT_ASSERT_EQUAL(_dbern->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_dbeta->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dbin->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dcat->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_dchisqr->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_ddirch->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_ddexp->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dexp->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_df->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dgamma->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dgengamma->npar(), 3U);
    CPPUNIT_ASSERT_EQUAL(_dhyper->npar(), 4U);
    CPPUNIT_ASSERT_EQUAL(_dinterval->npar(), 2U);    
    CPPUNIT_ASSERT_EQUAL(_dlnorm->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dlogis->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dmnorm->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dmnormvc->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dmt->npar(), 3U);
    CPPUNIT_ASSERT_EQUAL(_dmulti->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dchisqr->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_dnegbin->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dgampois->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dnorm->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dnt->npar(), 3U);
    CPPUNIT_ASSERT_EQUAL(_dpar->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dpois->npar(), 1U);
    CPPUNIT_ASSERT_EQUAL(_dround->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dsample->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dsum->npar(), 0U);
    CPPUNIT_ASSERT_EQUAL(_dt->npar(), 3U);
    CPPUNIT_ASSERT_EQUAL(_dunif->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dweib->npar(), 2U);
    CPPUNIT_ASSERT_EQUAL(_dwish->npar(), 2U);

}

void BugsDistTest::name()
{
    CPPUNIT_ASSERT_EQUAL(string("dbern"), _dbern->name());
    CPPUNIT_ASSERT_EQUAL(string("dbeta"), _dbeta->name());
    CPPUNIT_ASSERT_EQUAL(string("dbin"), _dbin->name());
    CPPUNIT_ASSERT_EQUAL(string("dcat"), _dcat->name());
    CPPUNIT_ASSERT_EQUAL(string("dchisqr"), _dchisqr->name());
    CPPUNIT_ASSERT_EQUAL(string("ddirch"), _ddirch->name());
    CPPUNIT_ASSERT_EQUAL(string("ddexp"), _ddexp->name());
    CPPUNIT_ASSERT_EQUAL(string("dexp"), _dexp->name());
    CPPUNIT_ASSERT_EQUAL(string("df"), _df->name());
    CPPUNIT_ASSERT_EQUAL(string("dgamma"), _dgamma->name());
    CPPUNIT_ASSERT_EQUAL(string("dgen.gamma"), _dgengamma->name());
    CPPUNIT_ASSERT_EQUAL(string("dhyper"), _dhyper->name());
    CPPUNIT_ASSERT_EQUAL(string("dinterval"), _dinterval->name());    
    CPPUNIT_ASSERT_EQUAL(string("dlnorm"), _dlnorm->name());
    CPPUNIT_ASSERT_EQUAL(string("dlogis"), _dlogis->name());
    CPPUNIT_ASSERT_EQUAL(string("dmnorm"), _dmnorm->name());
    CPPUNIT_ASSERT_EQUAL(string("dmnorm.vcov"), _dmnormvc->name());
    CPPUNIT_ASSERT_EQUAL(string("dmt"), _dmt->name());
    CPPUNIT_ASSERT_EQUAL(string("dmulti"), _dmulti->name());
    CPPUNIT_ASSERT_EQUAL(string("dnchisqr"), _dnchisqr->name());
    CPPUNIT_ASSERT_EQUAL(string("dnegbin"), _dnegbin->name());
    CPPUNIT_ASSERT_EQUAL(string("dgampois"), _dgampois->name());
    CPPUNIT_ASSERT_EQUAL(string("dnorm"), _dnorm->name());
    CPPUNIT_ASSERT_EQUAL(string("dnt"), _dnt->name());
    CPPUNIT_ASSERT_EQUAL(string("dpar"), _dpar->name());
    CPPUNIT_ASSERT_EQUAL(string("dpois"), _dpois->name());
    CPPUNIT_ASSERT_EQUAL(string("dround"), _dround->name());
    CPPUNIT_ASSERT_EQUAL(string("dsample"), _dsample->name());
    CPPUNIT_ASSERT_EQUAL(string("dsum"), _dsum->name());
    CPPUNIT_ASSERT_EQUAL(string("dt"), _dt->name());
    CPPUNIT_ASSERT_EQUAL(string("dunif"), _dunif->name());
    CPPUNIT_ASSERT_EQUAL(string("dweib"), _dweib->name());
    CPPUNIT_ASSERT_EQUAL(string("dwish"), _dwish->name());

}

void BugsDistTest::alias()
{
    CPPUNIT_ASSERT_EQUAL(string(""), _dbern->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dbeta->alias());
    CPPUNIT_ASSERT_EQUAL(string("dbinom"), _dbin->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dcat->alias());
    CPPUNIT_ASSERT_EQUAL(string("dchisq"), _dchisqr->alias());
    CPPUNIT_ASSERT_EQUAL(string("ddirich"), _ddirch->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _ddexp->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dexp->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _df->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dgamma->alias());
    CPPUNIT_ASSERT_EQUAL(string("dggamma"), _dgengamma->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dhyper->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dinterval->alias());    
    CPPUNIT_ASSERT_EQUAL(string(""), _dlnorm->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dlogis->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dmnorm->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dmnormvc->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dmt->alias());
    CPPUNIT_ASSERT_EQUAL(string("dmultinom"), _dmulti->alias());
    CPPUNIT_ASSERT_EQUAL(string("dnchisq"), _dnchisqr->alias());
    CPPUNIT_ASSERT_EQUAL(string("dnbinom"), _dnegbin->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dgampois->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dnorm->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dnt->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dpar->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dpois->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dround->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dsample->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dsum->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dt->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dunif->alias());
    CPPUNIT_ASSERT_EQUAL(string("dweibull"), _dweib->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dwish->alias());
}

void BugsDistTest::rscalar_rpq(RScalarDist const *dist, 
			       vector<double const *> const &par)
{
    /*
      Simultaneous test of r, p, and q functions for distributions
      inheriting from RScalarDist.
    */
    unsigned int nsim = 100;
    unsigned int nsim2 = 10;
    unsigned int nquant = 10;
    
    CPPUNIT_ASSERT_MESSAGE(dist->name(), checkNPar(dist, par.size()));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par));

    for (unsigned int s = 0; s < nsim; ++s) {
	//Generate random variable from distribution
	double y = dist->r(par, _rng);
	//Pass to distribution function and then to quantile function
	double p = dist->p(y, par, true, false);
	CPPUNIT_ASSERT_MESSAGE(dist->name(), p >= 0 && p <= 1);
	double z = dist->q(p, par, true, false);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z, tol*abs(y));
	//Now do the same on a log scale
	double logp = dist->p(y, par, true, true);
	CPPUNIT_ASSERT_MESSAGE(dist->name(), logp <= 0);
	z = dist->q(logp, par, true, true);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z, tol*abs(y));
	//Using upper tail
	p = dist->p(y, par, false, false);
	CPPUNIT_ASSERT_MESSAGE(dist->name(), p >= 0 && p <= 1);
	z = dist->q(p, par, false, false);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z, tol*abs(y));
	//Upper tail on log scale
	logp = dist->p(y, par, false, true);
	CPPUNIT_ASSERT_MESSAGE(dist->name(), logp <= 0);
	z = dist->q(logp, par, false, true);
	CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z, tol*abs(y));
    }

    for (unsigned int s = 0; s < nsim2; ++s) {
	
	// Test truncated sampling
	for (unsigned int i = 0; i <= nquant; ++i) {
	    double plim1 = static_cast<double>(i)/nquant;
	    double lim1 = dist->q(plim1, par, true, false);

	    if (i > 0) {
		// Test sampling from lower tail
		double y = dist->randomSample(par, 0, &lim1, _rng);
		CPPUNIT_ASSERT_MESSAGE(dist->name(), y <= lim1);
		double logpy = dist->p(y, par, true, true);
		double z = dist->q(logpy, par, true, true);
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z,
						     tol * abs(y));
	    }
	    if (i < nquant) {
		// Test sampling from upper tail
		double y = dist->randomSample(par, &lim1, 0, _rng);
		CPPUNIT_ASSERT_MESSAGE(dist->name(), y >= lim1);
		double logpy = dist->p(y, par, false, true);
		double z = dist->q(logpy, par, false, true);
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z,
						     tol * abs(y));
	    }
	    for (unsigned int j = i+1; j <= nquant; ++j) {
		// Test sampling between lower and upper limits
		double plim2 = static_cast<double>(j)/nquant;
		double lim2 = dist->q(plim2, par, true, false);

		double y = dist->randomSample(par, &lim1, &lim2, _rng);
		CPPUNIT_ASSERT_MESSAGE(dist->name(), y >= lim1 && y <= lim2);
		double py = dist->p(y, par, true, false);
		double z = dist->q(py, par, true, false);
		CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), y, z,
						     tol * abs(y));
	    }
	}
    }

    rscalar_trunclik(dist, par);
}

static void scalar_trunclik_discrete(RScalarDist const *dist,
				     vector<double const *> const &par,
				     double bound, bool lower)
{
    // Test normalization of truncated likelihood for discrete variables
    // by enumeration
    
    CPPUNIT_ASSERT(dist->discrete());

    double delta = lower ? -1 : 1;
    double const *ll = 0;
    double const *uu = 0;
    if (lower) {
	uu = &bound;
    }
    else {
	ll = &bound;
    }
    double lik = 0;
    double likmax = dist->logDensity(bound, jags::PDF_FULL, par, ll, uu);

    double x = bound;
    for (unsigned int i = 0; i < 1000; ++i) {
	double y = dist->logDensity(x, jags::PDF_FULL, par,  ll, uu);
	lik += exp(y);
	if (y > likmax) {
	    likmax = y;
	}
	else if (y < likmax - 10) {
	    break;
	}
	x += delta;
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), 1.0, lik, 1e-2);
}

static void scalar_trunclik_cont(RScalarDist const *dist,
				 vector<double const *> const &par,
				 double bound, bool lower)
{
    // Test normalization of truncated likelihood for continuous
    // variables by approximate integration using the trapezoid method
    
    CPPUNIT_ASSERT(!dist->discrete());
    CPPUNIT_ASSERT(jags_finite(bound));
    
    //Ensure that the density is finite at the other boundary (ob)
    //If not then just return as this method will not work
    double ob = lower ? dist->l(par) : dist->u(par);
    if (jags_finite(ob)) {
	if (!jags_finite(dist->d(ob, jags::PDF_FULL, par, false))) {
	    return;
	}
    }

    //Create grid of points for evaluating the likelihood
    //If lower==false then the grid is in reverse order
    vector<double> x;
    unsigned int N = 1000;
    if (jags_finite(ob)) {
	x.push_back(ob);
    }
    double pmax = dist->p(bound, par, lower, false);
    for (unsigned int i = 1; i < N; ++i) {
	x.push_back(dist->q(i*pmax/N, par, lower, false));
    }
    x.push_back(bound);

    // Calculate likelihood by trapezoid method
    // FIXME: WE could probably do better
    double lik = 0;
    double const *ll = lower ? 0 : &bound;
    double const *uu = lower ? &bound : 0;
    double y0 = exp(dist->logDensity(x[0], jags::PDF_FULL, par, ll, uu)); 
    for (unsigned int i = 1; i < x.size(); ++i) {
	double y = exp(dist->logDensity(x[i], jags::PDF_FULL, par, ll, uu));
	lik += (x[i] - x[i-1]) * (y + y0) / 2;
	y0 = y;
    }
    if (!lower) {
	lik = -lik; //adjust for reverse grid
    }
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), 1.0, lik, 1.1e-2);
}


void BugsDistTest::rscalar_trunclik(RScalarDist const *dist, 
				    vector<double const *> const &par)
{
    /*
      Test likelihood calculations for truncated distributions
    */
    CPPUNIT_ASSERT_MESSAGE(dist->name(), checkNPar(dist, par.size()));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par));
    
    double ll = dist->q(0.1, par, true, false);
    double uu = dist->q(0.9, par, true, false);

    if (dist->discrete()) {
	scalar_trunclik_discrete(dist, par, ll, true);
	scalar_trunclik_discrete(dist, par, ll, false);
	scalar_trunclik_discrete(dist, par, uu, true);
	scalar_trunclik_discrete(dist, par, uu, false);
    }
    else {
	scalar_trunclik_cont(dist, par, ll, true);
	scalar_trunclik_cont(dist, par, ll, false);
	scalar_trunclik_cont(dist, par, uu, true);
	scalar_trunclik_cont(dist, par, uu, false);
    }
}


static vector<double const *> mkPar(double const &v1)
{
    vector<double const *> par;
    par.push_back(&v1);
    return par;
}

static vector<double const *> mkPar(double const &v1, double const &v2)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    return par;
}

static vector<double const *> mkPar(double const &v1, double const &v2,
				    double const &v3)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    par.push_back(&v3);
    return par;
}

static vector<double const *> mkPar(double const &v1, double const &v2,
				    double const &v3, double const &v4)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    par.push_back(&v3);
    par.push_back(&v4);
    return par;
}

void BugsDistTest::rscalar()
{
    rscalar_rpq(_dbeta, mkPar(0.3, 0.5));
    rscalar_rpq(_dbeta, mkPar(2, 6));
    
    rscalar_rpq(_dbin, mkPar(0.1, 10));
    rscalar_rpq(_dbin, mkPar(0.9, 1));

    rscalar_rpq(_dchisqr, mkPar(1));
    rscalar_rpq(_dchisqr, mkPar(2));
    rscalar_rpq(_dchisqr, mkPar(3));
    rscalar_rpq(_dchisqr, mkPar(10));

    rscalar_rpq(_ddexp, mkPar(2, 0.1));
    rscalar_rpq(_ddexp, mkPar(-5, 2.8));

    
    rscalar_rpq(_dexp, mkPar(0.1));
    rscalar_rpq(_dexp, mkPar(7));

    rscalar_rpq(_df, mkPar(1, 1));
    rscalar_rpq(_df, mkPar(3, 7));

    rscalar_rpq(_dgamma, mkPar(1, 0.1));
    rscalar_rpq(_dgamma, mkPar(0.1, 3));

    rscalar_rpq(_dgengamma, mkPar(1, 0.1, 3));
    rscalar_rpq(_dgengamma, mkPar(0.1, 3, 0.5));

    for (int i = 1; i <= 10; ++i) {
	for (int j = 0; j <= i; ++j ) {
	    rscalar_rpq(_dhyper, mkPar(i, i, j, 1));
	    rscalar_rpq(_dhyper, mkPar(i, i, j, 0.5));
	    rscalar_rpq(_dhyper, mkPar(i, i, j, 3));
	}
    }
    rscalar_rpq(_dhyper, mkPar(10, 8, 7, 0.5));

    rscalar_rpq(_dlnorm, mkPar(-0.5, 33));
    rscalar_rpq(_dlogis, mkPar(-7, 0.3));
    rscalar_rpq(_dnchisqr, mkPar(5, 1.5));
    rscalar_rpq(_dnorm, mkPar(3, 10));
    rscalar_rpq(_dnt, mkPar(0, 1, 4));
    rscalar_rpq(_dnt, mkPar(-5, 3, 1));
    rscalar_rpq(_dpar, mkPar(1.3, 1));

    rscalar_rpq(_dpois, mkPar(0.12));
    rscalar_rpq(_dpois, mkPar(10.7));
    rscalar_rpq(_dt, mkPar(18, 1.7, 1));
    rscalar_rpq(_dt, mkPar(0, 1, 4));
    rscalar_rpq(_dweib, mkPar(3, 10));
    rscalar_rpq(_dweib, mkPar(0.9, 0.9));
}

void BugsDistTest::kl_scalar(ScalarDist const *dist, 
			     vector<double const *> const &par0,
			     vector<double const *> const &par1)
{
    /* 
       Test closed form expression for Kullback-Leibler divergence
       against Monte Carlo estimate
    */
    
    //Sanity checks
    CPPUNIT_ASSERT_EQUAL_MESSAGE(dist->name(), par0.size(), par1.size());
    CPPUNIT_ASSERT_MESSAGE(dist->name(), checkNPar(dist, par0.size()));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par0));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par1));
    
    unsigned int nrep = 1000L;

    //Monte Carlo estimate of Kullback-Leibler divergence
    
    vector<double> y(nrep);
    double ysum = 0;
    for (unsigned int i = 0; i < nrep; ++i) {
	y[i] = dist->KL(par0, par1, 0, 0, _rng, 100);
	ysum += y[i];
    }
    double ymean = ysum/nrep;
    double ysum2 = 0;
    for (unsigned int i = 0; i < nrep; ++i) {
	ysum2 += (y[i] - ymean) * (y[i] - ymean);
    }
    double ysd = sqrt(ysum2/(nrep - 1)); //Standard deviation
    //Standard error of the mean
    double ysem = ysd/sqrt(static_cast<double>(nrep)); 
    
    //Closed form expression for Kullback-Leibler divergence
    double z = dist->KL(par0, par1);
    CPPUNIT_ASSERT(z != JAGS_NA);

    //Difference between closed-form and Monte Carlo estimates
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(dist->name(), z, ymean,
					 max(5*ysem, tol));
}

void BugsDistTest::kl()
{
    /* Kullback-Leibler test: see kl_scalar for details */
    
    kl_scalar(_dbern, mkPar(0.6), mkPar(0.5));
    kl_scalar(_dbern, mkPar(0.1), mkPar(0.9));

    kl_scalar(_dbin, mkPar(0.5, 10), mkPar(0.9, 10));
    kl_scalar(_dbin, mkPar(0.1, 5), mkPar(0.2, 5));

    kl_scalar(_dbeta, mkPar(2, 3), mkPar(2.5, 2.8));
    kl_scalar(_dbeta, mkPar(0.5, 1), mkPar(0.7, 1));

    kl_scalar(_dchisqr, mkPar(1), mkPar(1.5));
    kl_scalar(_dchisqr, mkPar(5), mkPar(7));

    kl_scalar(_ddexp, mkPar(0, 1), mkPar(-1, 3));
    kl_scalar(_ddexp, mkPar(5, 0.2), mkPar(7, 0.25));

    kl_scalar(_dexp, mkPar(10), mkPar(12));
    kl_scalar(_dexp, mkPar(0.5), mkPar(0.7));

    kl_scalar(_dgamma, mkPar(3,2), mkPar(3.5, 1.8));
    kl_scalar(_dgamma, mkPar(1, 0.5), mkPar(2, 0.1));

    kl_scalar(_dgengamma, mkPar(2,5,3), mkPar(2.5, 4.8, 2.8));
    kl_scalar(_dgengamma, mkPar(1, 3, 0.5), mkPar(1.1, 2.2, 0.7));

    kl_scalar(_dhyper, mkPar(5,5,6,0.2), mkPar(5,5,6,0.4));
    kl_scalar(_dhyper, mkPar(20, 10, 12, 0.8), mkPar(20,10, 12, 0.6));
    
    kl_scalar(_dlnorm, mkPar(0, 1), mkPar(0.5, 1.5));
    kl_scalar(_dlnorm, mkPar(3, 2), mkPar(5, 2.8));

    //Can only calculate for negative binomial when sizes are equal
    kl_scalar(_dnegbin, mkPar(0.3, 10), mkPar(0.6, 10));
    kl_scalar(_dnegbin, mkPar(0.1, 3), mkPar(0.2, 3));
    
    kl_scalar(_dnorm, mkPar(0, 1), mkPar(-3, 2));
    kl_scalar(_dnorm, mkPar(2, 0.3), mkPar(2.5, 0.1));

    kl_scalar(_dpois, mkPar(3), mkPar(5));
    kl_scalar(_dpois, mkPar(12), mkPar(10));

    kl_scalar(_dunif, mkPar(-5,10), mkPar(-8, 12));
    kl_scalar(_dunif, mkPar(2,3), mkPar(1,5));

    /* FIXME
       Need tests for DCat, DDirch
       Should implemnt DMNorm::KL
    */
}    

static double superror(RScalarDist const *dist, unsigned int N,
		       vector<double const *> const &par, jags::RNG *rng)
{
    //Calculate test statistic for DKW test
    
    multiset<double> xmset;
    
    for (unsigned int i =0; i < N; ++i) {
	double x = dist->r(par, rng);
	x = fprec(x, 12);
	xmset.insert(x);
    }
    double fhat=0, delta=0;
    for (multiset<double>::const_iterator q = xmset.begin();
	 q != xmset.end(); q = xmset.upper_bound(*q))
    {
	double f = dist->p(*q, par, true, false);
	fhat += static_cast<double>(xmset.count(*q))/N;
	delta = max(delta, abs(fhat - f));
    }
    return delta;
}

static double pdkwbound(double n, double t) {
    return 2*exp(-2*n*t*t);
}

static double qdkwbound(double n, double p) {
    return sqrt(log(p/2)/(-2*n));
}

void BugsDistTest::dkwtest(RScalarDist const *dist,
			   vector<double const *> const &par,
			   unsigned int N, double pthresh)
{
    /*
      Test using the Dvoretzky-Kiefer-Wolfowitz (1956) bound on the
      difference between the empirical distribution function and the
      theoretical one, with a tight constant derived by Massart (1990)
    */
    
    //Sanity checks
    CPPUNIT_ASSERT_MESSAGE(dist->name(), checkNPar(dist, par.size()));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par));

    double s = superror(dist, N, par, _rng);
    string msg = dist->name();
    if (s >= qdkwbound(N, pthresh)) {
	ostringstream ostr;
	ostr << "Distribution " << dist->name() << "(";
	for (unsigned int i = 0; i < par.size(); ++i) {
	    if (i > 0) ostr << ", ";
	    ostr << *par[i];
	}
	ostr << ")\n";
	ostr << "supremum error = " << fprec(s, 2)
	     << " with p-value="
	     << min(1.0, fround(pdkwbound(N,s), 4));
	CPPUNIT_FAIL(ostr.str());
    }

	
}

//FIXME: We should have a atest for simulation from a truncated distribution

void BugsDistTest::dkw()
{
    //DKW test of random number generation and distribution function
    //See dkwtest for details
    
    //CPPUNIT_ASSERT_MESSAGE("dkw check", false);

    dkwtest(_dbeta, mkPar(1,1));
    dkwtest(_dbeta, mkPar(2,1));
    dkwtest(_dbeta, mkPar(1,2));
    dkwtest(_dbeta, mkPar(2,2));
    dkwtest(_dbeta, mkPar(0.2,0.2));

    dkwtest(_dbin, mkPar(0.2, 1));
    dkwtest(_dbin, mkPar(0.2, 2));
    dkwtest(_dbin, mkPar(0.2, 100));
    dkwtest(_dbin, mkPar(0.2, 1e4));
    dkwtest(_dbin, mkPar(0.8, 1));
    dkwtest(_dbin, mkPar(0.8, 100));
    dkwtest(_dbin, mkPar(0.999, 100));

    dkwtest(_dchisqr, mkPar(1));
    dkwtest(_dchisqr, mkPar(10));

    dkwtest(_ddexp, mkPar(0, 1));
    dkwtest(_ddexp, mkPar(0, 2));
    dkwtest(_ddexp, mkPar(-3, 0.5));

    dkwtest(_dexp, mkPar(1));
    dkwtest(_dexp, mkPar(10));
    dkwtest(_dexp, mkPar(0.01));
    
    dkwtest(_df, mkPar(1, 1));
    dkwtest(_df, mkPar(1, 10));
    dkwtest(_df, mkPar(10, 10));
    dkwtest(_df, mkPar(30, 3));

    dkwtest(_dgamma, mkPar(0.1, 1));
    dkwtest(_dgamma, mkPar(0.2, 1));
    dkwtest(_dgamma, mkPar(10, 1));
    dkwtest(_dgamma, mkPar(20, 1));
    
    dkwtest(_dgengamma, mkPar(1, 0.1, 0.5));
    dkwtest(_dgengamma, mkPar(0.1, 7, 2));

    dkwtest(_dpois, mkPar(0.095));
    dkwtest(_dpois, mkPar(0.95));
    dkwtest(_dpois, mkPar(9.5));
    dkwtest(_dpois, mkPar(95));

    dkwtest(_dnegbin, mkPar(0.2, 1));
    dkwtest(_dnegbin, mkPar(0.2, 2));
    dkwtest(_dnegbin, mkPar(0.2, 100));
    dkwtest(_dnegbin, mkPar(0.2, 1e4));
    dkwtest(_dnegbin, mkPar(0.8, 1));
    dkwtest(_dnegbin, mkPar(0.8, 100));
    dkwtest(_dnegbin, mkPar(0.999, 100));

    dkwtest(_dnorm, mkPar(0, 1));
    dkwtest(_dnorm, mkPar(5, 3));

    //Central hypergeometric
    dkwtest(_dhyper, mkPar(40, 30, 10, 1));
    dkwtest(_dhyper, mkPar(40, 3, 20, 1));
    dkwtest(_dhyper, mkPar(40, 3, 20, 1));
    dkwtest(_dhyper, mkPar(6, 3, 2, 1));
    dkwtest(_dhyper, mkPar(5, 3, 2, 1));
    dkwtest(_dhyper, mkPar(4, 3, 2, 1));

    //Large sample test for central hypergeometric
    dkwtest(_dhyper, mkPar(60,100,50, 1));
    dkwtest(_dhyper, mkPar(600, 1000, 500, 1));

    //Non-central hypergeometric
    dkwtest(_dhyper, mkPar(40, 30, 10, 2));
    dkwtest(_dhyper, mkPar(40, 3, 20, 0.5));
    dkwtest(_dhyper, mkPar(40, 3, 20, 0.3));
    dkwtest(_dhyper, mkPar(6, 3, 2, 5));
    dkwtest(_dhyper, mkPar(5, 3, 2, 2));
    dkwtest(_dhyper, mkPar(4, 3, 2, 2.5));

    dkwtest(_dlnorm, mkPar(0, 1));
    dkwtest(_dlnorm, mkPar(3, 3));
    dkwtest(_dlnorm, mkPar(-8, 0.2));
    
    dkwtest(_dlogis, mkPar(0, 1));
    dkwtest(_dlogis, mkPar(4, 0.5));

    dkwtest(_dnt, mkPar(10, 2, 3));
    dkwtest(_dnt, mkPar(-7, 0.5, 1.5));

    dkwtest(_dpar, mkPar(5,2));
    
    dkwtest(_dt, mkPar(0,1,1));
    dkwtest(_dt, mkPar(-3, 0.4, 1));
    dkwtest(_dt, mkPar(0, 2, 10));
    dkwtest(_dt, mkPar(5, 1, 40));

    dkwtest(_dweib, mkPar(1, 1));
    dkwtest(_dweib, mkPar(2, 2));
    dkwtest(_dweib, mkPar(0.3, 0.5));
}
    
