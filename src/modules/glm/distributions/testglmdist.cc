#include <config.h>

#include "testglmdist.h"

#include "DScaledGamma.h"
#include "DScaledWishart.h"
#include "DOrderedLogit.h"

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

#define F77_DPOTRF F77_FUNC(dpotrf,DPOTRF)
#define F77_DPOTRI F77_FUNC(dpotri, DPOTRI)

using jags::RScalarDist;

extern "C" {
    void F77_DPOTRF (const char *uplo, const int *n, double *a,
		     const int *lda, const int *info);
    void F77_DPOTRI (const char *uplo, const int *n, double *a,
		     const int *lda, const int *info);

}

void GLMDistTest::setUp() {

    _rng = new jags::base::MersenneTwisterRNG(1234567, 
					      jags::KINDERMAN_RAMAGE);

    _dscaled_gamma = new jags::glm::DScaledGamma();
    _dscaled_wishart = new jags::glm::DScaledWishart();
    _dordered_logit = new jags::glm::DOrderedLogit();
}

void GLMDistTest::tearDown() {

    delete _rng;

    delete _dscaled_gamma;
    delete _dscaled_wishart;
    delete _dordered_logit;
}

void GLMDistTest::npar()
{
    //CPPUNIT_ASSERT_MESSAGE("npar check", false);

    CPPUNIT_ASSERT_EQUAL(_dscaled_gamma->npar(), 2UL);
    CPPUNIT_ASSERT_EQUAL(_dscaled_wishart->npar(), 2UL);
    CPPUNIT_ASSERT_EQUAL(_dordered_logit->npar(), 2UL);
}

void GLMDistTest::name()
{
    CPPUNIT_ASSERT_EQUAL(string("dscaled.gamma"), _dscaled_gamma->name());
    CPPUNIT_ASSERT_EQUAL(string("dscaled.wishart"), _dscaled_wishart->name());
    CPPUNIT_ASSERT_EQUAL(string("dordered.logit"), _dordered_logit->name());
}

void GLMDistTest::alias()
{
    CPPUNIT_ASSERT_EQUAL(string(""), _dscaled_gamma->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dscaled_wishart->alias());
    CPPUNIT_ASSERT_EQUAL(string(""), _dordered_logit->alias());
}

void GLMDistTest::rscalar_rpq(RScalarDist const *dist, 
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


void GLMDistTest::rscalar_trunclik(RScalarDist const *dist, 
				    vector<double const *> const &par)
{
    /*
      Test likelihood calculations for truncated distributions
    */
    CPPUNIT_ASSERT_MESSAGE(dist->name(), checkNPar(dist, par.size()));
    CPPUNIT_ASSERT_MESSAGE(dist->name(), dist->checkParameterValue(par));
    
    double ll = dist->q(0.1, par, true, false);
    double uu = dist->q(0.9, par, true, false);

    scalar_trunclik_cont(dist, par, ll, true);
    scalar_trunclik_cont(dist, par, ll, false);
    scalar_trunclik_cont(dist, par, uu, true);
    scalar_trunclik_cont(dist, par, uu, false);
}


/*
static vector<double const *> mkPar(double const &v1)
{
    vector<double const *> par;
    par.push_back(&v1);
    return par;
}
*/

static vector<double const *> mkPar(double const &v1, double const &v2)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    return par;
}

/*
static vector<double const *> mkPar(double const &v1, double const &v2,
				    double const &v3)
{
    vector<double const *> par;
    par.push_back(&v1);
    par.push_back(&v2);
    par.push_back(&v3);
    return par;
}
*/

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

void GLMDistTest::rscalar()
{
    rscalar_rpq(_dscaled_gamma, mkPar(0.3, 1));
    rscalar_rpq(_dscaled_gamma, mkPar(5.0, 2));
    rscalar_rpq(_dscaled_gamma, mkPar(2.0, 3));
    rscalar_rpq(_dscaled_gamma, mkPar(0.7, 4));
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

void GLMDistTest::dkwtest(RScalarDist const *dist,
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

void GLMDistTest::dkw()
{
    //DKW test of random number generation and distribution function
    //See dkwtest for details
    
    //CPPUNIT_ASSERT_MESSAGE("dkw check", false);

    dkwtest(_dscaled_gamma, mkPar(1,1));
    dkwtest(_dscaled_gamma, mkPar(1,2));
    dkwtest(_dscaled_gamma, mkPar(1,3));
    dkwtest(_dscaled_gamma, mkPar(1,4));

    dkwtest(_dscaled_gamma, mkPar(3,1));
    dkwtest(_dscaled_gamma, mkPar(4,2));
    dkwtest(_dscaled_gamma, mkPar(5,3));
    dkwtest(_dscaled_gamma, mkPar(6,4));
}


    
static bool inverse_spd (double *A, int n)
{
    /* invert n x n symmetric positive definite matrix A*/

    int info = 0;
    F77_DPOTRF ("L", &n, A, &n, &info);
    if (info < 0) {
	CPPUNIT_FAIL("Illegal argument in inverse_spd");
    }
    else if (info > 0) {
	CPPUNIT_FAIL("Cannot invert matrix: not positive definite");
    }
    F77_DPOTRI ("L", &n, A, &n, &info); 

    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < i; ++j) {
	    A[i*n + j] = A[j*n + i];
	}
    }

    if (info != 0) {
	CPPUNIT_FAIL("Unable to invert symmetric positive definite matrix");
    }
    return true;
}

void GLMDistTest::dkw_swish(vector<double> const &S,
			    double df, unsigned int r, unsigned int c,
			    unsigned int N, double pthresh)
{
    //Special version of the DKW test for the scaled Wishart distribution
    //We invert the matrix and look at individual elements of the inverse:
    //Diagonal elements should be scaled t
    //Off-diagonal elements should be scaled beta
    
    multiset<double> xmset;

    unsigned int ndim = S.size();
    unsigned int length = ndim * ndim;
    vector<double> x(length);

    //Set up parameters
    vector<double const*> par;
    par.push_back(&S[0]);
    par.push_back(&df);

    //Set up dimensions
    vector<vector<unsigned long> > dims;
    dims.push_back(vector<unsigned long>(1, ndim));
    dims.push_back(vector<unsigned long>(1, 1));

    //Sanity checks
    CPPUNIT_ASSERT_MESSAGE("scaled.wishart",
			   checkNPar(_dscaled_wishart, par.size()));
    CPPUNIT_ASSERT_MESSAGE("scaled.wishart",
			   _dscaled_wishart->checkParameterDim(dims));
    CPPUNIT_ASSERT_MESSAGE("scaled.wishart",
			   _dscaled_wishart->checkParameterValue(par, dims));
			   
    //Calculate maximum error for DKW test.
    for (unsigned int i = 0; i < N; ++i) {
	_dscaled_wishart->randomSample(&x[0], par, dims, _rng);

	inverse_spd(&x[0], ndim);
	
	double y = x[r*ndim + c];
	if (r == c) {
	    //Standard deviation
	    y = sqrt(y) / S[r]; 
	}
	else {
	    //Correlation
	    y /= sqrt(x[r*ndim+r]*x[c*ndim+c]);
	    y = (y+1)/2;
	}
	y = fprec(y, 12);
	xmset.insert(y);
    }
    double fhat=0, delta=0;
    for (multiset<double>::const_iterator q = xmset.begin();
	 q != xmset.end(); q = xmset.upper_bound(*q))
    {
	double f;
	if (r == c) {
	    //For diagonal elements, use half-t distribution
	    f = 2 * pt(*q, df, true, false) - 1;
	}
	else {
	    //For off-diagonal elements, use beta-distribution
	    f = pbeta(*q, df/2, df/2, true, false);
	}
	fhat += static_cast<double>(xmset.count(*q))/N;
	delta = max(delta, abs(fhat - f));
    }

    if (delta >= qdkwbound(N, pthresh)) {
	ostringstream ostr;
	ostr << "Distribution scaled.wishart(scale=c(";
	for (unsigned int i = 0; i < S.size(); ++i) {
	    if (i > 0) ostr << ", ";
	    ostr << S[i];
	}
	ostr << "), df=" << df << ")\n";
	ostr << "Row = " << r << ", Column = " << c << "\n";
	ostr << "supremum error = " << fprec(delta, 2)
	     << " with p-value="
	     << min(1.0, fround(pdkwbound(N,delta), 4));
	CPPUNIT_FAIL(ostr.str());
    }
}

void GLMDistTest::scaled_wishart() {

    /* Checks that the elements of the inverse matrix have the correct
       marginal distribution */
    
    //2x2 matrix, unscaled
    
    vector<double> S(2,1);

    for (unsigned int r = 0; r < S.size(); ++r) {
	for (unsigned int c = 0; c < S.size(); ++c) {
	    dkw_swish(S, 2, r, c);
	    dkw_swish(S, 1, r, c);
	    dkw_swish(S, 4.8, r, c);
	    dkw_swish(S, 10, r, c);
	}
    }

    //3x3 matrix, scaled
    
    S.clear();
    S.push_back(2);
    S.push_back(1.3);
    S.push_back(0.5);
    
    for (unsigned int r = 0; r < S.size(); ++r) {
	for (unsigned int c = 0; c < S.size(); ++c) {
	    dkw_swish(S, 2, r, c);
	    dkw_swish(S, 1, r, c);
	    dkw_swish(S, 1.3, r, c);
	    dkw_swish(S, 3.7, r, c);
	    dkw_swish(S, 5.5, r, c);
	}
    }
	

}
