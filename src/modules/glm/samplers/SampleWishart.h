#ifndef SAMPLE_WISHART_H_
#define SAMPLE_WISHART_H_

namespace jags {

    class RNG;

    namespace glm {
	/** @short Sample from Wishart distribution */
	void sampleWishart(double *X, int length, double const *R,
			   double df, int nrow, RNG *rng);
    }
}

#endif /* SAMPLE_WISHART_H_ */
