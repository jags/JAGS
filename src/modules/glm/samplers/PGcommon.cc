#include "PG.h"

#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <JRmath.h>

using jags::RNG;

static double rigauss_body(double imu, double lambda, double t, RNG *rng)
{
    // Sample truncated IG(mu, lambda) I(0,t) using accept-reject sampling
    // with a truncated inverse chi-square as proposal. This is efficient
    // when the truncation point t is in the body of the distribution.
    // Here mu = 1/imu
    
    // Rescale the problem to sample Z ~ IG(mu/lambda, 1) I(0, t/lambda)
    imu *= lambda;
    t /= lambda;
    
    double alpha, Z;
    do {
	Z = lnormal(1/sqrt(t), rng); 
	Z = 1/(Z*Z); // Z ~ 1/dchisq(1) truncated to (0,t)
	alpha = exp(- imu * imu * Z / 2.0); //acceptance probability
    }
    while (rng->uniform() >= alpha);
    
    return Z * lambda; // Rescale the problem to sample X ~ IG(mu, lambda) I(0,t)
}

static double rigauss_tail(double mu, double lambda, double t, RNG *rng)
{
    // Generate truncated IG(mu, lambda) I(0, t) by rejection sampling
    // This is efficient when the truncation point t is in the tail of
    // the distribution.
    
    double X;
    do {
	// Repeatedly generate X ~ IG(mu, lambda) according to
	// Devroye (1986) until X < t
	double Y = rng->normal();
	Y *= Y;
	double muY = mu * Y;
	X = mu + mu * (muY - sqrt((4 * lambda + muY) * muY)) / (2 * lambda);

	if (rng->uniform() > mu / (mu + X)) {
	    X = mu * mu / X;
	}
    }
    while (X >= t);

    return X;
}

namespace jags {
    namespace glm {
	
	double rigauss(double imu, double lambda, double t, RNG *rng)
	{
	    // Sample from the inverse Gaussian distribution with
	    // mean mu=1/imu and shape lambda, truncated to the interval (0, t)

	    if (imu * t < 1.0) {
		return rigauss_body(imu, lambda, t, rng);
	    }
	    else {
		return rigauss_tail(1/imu, lambda, t, rng);
	    }
	}
	
    }
}

