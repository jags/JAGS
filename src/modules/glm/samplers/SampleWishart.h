#ifndef SAMPLE_WISHART_H_
#define SAMPLE_WISHART_H_

namespace jags {

    struct RNG;

    namespace glm {
	/** @short Sample from Wishart distribution */
	void sampleWishart(double *X, unsigned long length, double const *R,
			   double df, unsigned long nrow, RNG *rng);
    }
}

#endif /* SAMPLE_WISHART_H_ */
