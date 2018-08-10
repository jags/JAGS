#ifndef PG_H_
#define PG_H_

namespace jags {

    struct RNG;
    
    namespace glm {
	/**
	 * Sample from a inverse Gaussian distribution truncated to (0,t)
	 *
	 * @param imu Reciprocal of the mean parameter
	 * @param lambda Shape parameter
	 * @param t Truncation point
	 * @param rng Random number generator
	 */
	double rigauss(double imu, double lambda, double t, RNG *rng);
    }
}

#endif /* PG_H_ */
