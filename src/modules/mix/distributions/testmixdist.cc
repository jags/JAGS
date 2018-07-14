#include "testmixdist.h"

#include "DBetaBin.h"
#include "DNormMix.h"

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
using jags::PDF_FULL;

void MixDistTest::setUp() {

    _rng = new jags::base::MersenneTwisterRNG(1234567, 
					      jags::KINDERMAN_RAMAGE);

    _dbetabin = new jags::mix::DBetaBin();
    _dnormmix = new jags::mix::DNormMix();
}

void MixDistTest::tearDown() {

    delete _rng;

    delete _dbetabin;
    delete _dnormmix;
}

void MixDistTest::npar()
{
    //CPPUNIT_ASSERT_MESSAGE("npar check", false);

    CPPUNIT_ASSERT_EQUAL(_dbetabin->npar(), 3UL);
    CPPUNIT_ASSERT_EQUAL(_dnormmix->npar(), 3UL);
}

void MixDistTest::name()
{
    CPPUNIT_ASSERT_EQUAL(string("dbetabin"), _dbetabin->name());
    CPPUNIT_ASSERT_EQUAL(string("dnormmix"), _dnormmix->name());
}

void MixDistTest::alias()
{
    CPPUNIT_ASSERT_EQUAL(string("dbetabinom"), _dbetabin->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dnormmix->alias());
}

void MixDistTest::rscalar_rpq(RScalarDist const *dist, 
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

    //rscalar_trunclik(dist, par);
}

/*
static vector<double const *> mkPar(double const &v1)
{
    vector<double const *> par;
    par.push_back(&v1);
    return par;
}
*/

/*
static vector<double const *> mkPar(double const &v1, double const &v2)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    return par;
}
*/

static vector<double const *> mkPar(double const &v1, double const &v2,
				    double const &v3)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    par.push_back(&v3);
    return par;
}

/*
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
*/

void MixDistTest::rscalar()
{
    rscalar_rpq(_dbetabin, mkPar(1, 1, 5));
    rscalar_rpq(_dbetabin, mkPar(2, 1.3, 8));
    rscalar_rpq(_dbetabin, mkPar(1.8, 5, 8));
    rscalar_rpq(_dbetabin, mkPar(1, 1, 20));
    rscalar_rpq(_dbetabin, mkPar(2, 2, 30));
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

void MixDistTest::dkwtest(RScalarDist const *dist,
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

void MixDistTest::dkw()
{
    //DKW test of random number generation and distribution function
    //See dkwtest for details
    
    //CPPUNIT_ASSERT_MESSAGE("dkw check", false);

    dkwtest(_dbetabin, mkPar(1,1, 4));
}

static double normmix_mean(vector<double const *> const &par,
			   vector<unsigned long> const &lengths)
{
    double const *mu = par[0];
    double const *pi = par[2];
    
    double sumpi = 0, mean = 0;
    for (unsigned int i = 0; i < lengths[0]; ++i) {
	sumpi += pi[i];
	mean += pi[i] * mu[i];
    }
    return mean/sumpi;
}


static double normmix_var(vector<double const *> const &par,
			  vector<unsigned long> const &len)
{
    double const *mu = par[0];
    double const *tau = par[1];
    double const *pi = par[2];
    
    double sumpi = 0, mean = 0;
    for (unsigned int i = 0; i < len[0]; ++i) {
	sumpi += pi[i];
	mean += pi[i] * mu[i];
    }
    mean /= sumpi;
    
    double var = 0;
    for (unsigned int i = 0; i < len[0]; ++i) {
	var += pi[i] * ((mu[i] - mean)*(mu[i] - mean) + 1/tau[i]);
    }
    var /= sumpi;
    
    return var;
}

void MixDistTest::test_mean_normmix(vector<double const *> const &par,
				    vector<unsigned long> const &len,
				    unsigned int N)
{
    CPPUNIT_ASSERT_MESSAGE(_dnormmix->name(),
			   _dnormmix->checkParameterLength(len));
    CPPUNIT_ASSERT_MESSAGE(_dnormmix->name(),
			   _dnormmix->checkParameterValue(par, len));
    
    double x;
    double xmean = 0;
    for (unsigned int i = 0; i < N; ++i) {
	_dnormmix->randomSample(&x, par, len, 0, 0, _rng);
	xmean += x;
    }
    xmean /= N;

    double zmean = normmix_mean(par, len);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(_dnormmix->name(), zmean, xmean,
					 0.01);
}

void MixDistTest::test_var_normmix(vector<double const *> const &par,
				   vector<unsigned long> const &len,
				   unsigned int N)
{
    CPPUNIT_ASSERT_MESSAGE(_dnormmix->name(),
			   _dnormmix->checkParameterLength(len));
    CPPUNIT_ASSERT_MESSAGE(_dnormmix->name(),
			   _dnormmix->checkParameterValue(par, len));
    
    double x;
    double xmean = 0;
    for (unsigned int i = 0; i < N; ++i) {
	_dnormmix->randomSample(&x, par, len, 0, 0, _rng);
	xmean += x;
    }
    xmean /= N;

    double xvar = 0;
    for (unsigned int i = 0; i < N; ++i) {
	_dnormmix->randomSample(&x, par, len, 0, 0, _rng);
	xvar += (x - xmean) * (x - xmean);
    }
    xvar /= N;

    double zvar = normmix_var(par, len);
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(_dnormmix->name(),
					 zvar, xvar, 0.1);
}

void MixDistTest::normmix()
{
    unsigned int N = 100000;
    
    //Construct parameter vector of length 2
    vector<double> mu(2), tau(2), pi(2);
    vector<double const *> par2;
    par2.push_back(&mu[0]);
    par2.push_back(&tau[0]);
    par2.push_back(&pi[0]);
    
    //Construct length vector corresponding to par2
    vector<unsigned long> len2(3, 2);

    //Set up example
    mu[0] = -1; mu[1] = 3;
    tau[0] = 0.1; tau[1] = 0.2;
    pi[0] = 8; pi[1] = 2;
    test_mean_normmix(par2, len2, N);
    test_var_normmix(par2, len2, N);
    
    mu[0] = 1; mu[1] = 1;
    tau[0] = 1; tau[1] = 1;
    test_mean_normmix(par2, len2, N);
    test_var_normmix(par2, len2, N);
}
