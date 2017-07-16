#ifndef RE_METHOD2_H_
#define RE_METHOD2_H_

#include <sampler/MutableSampleMethod.h>

extern "C" {
#include <cholmod.h>
}

#include <vector>

namespace jags {

    struct RNG;
    class GraphView;
    class SingletonGraphView;

namespace glm {

    class GLMMethod;
    class Outcome;

    /**
     * @short Sample method for random effects and their variance.
     */
    class REMethod2 : public MutableSampleMethod {
      protected:
	SingletonGraphView const *_tau;
	GraphView const *_eps;
	std::vector<Outcome*> const &_outcomes;
	cholmod_sparse const *_x;
	const unsigned int _chain;
	cholmod_dense *_z;
	std::vector<unsigned int> _indices;
      public:
	/**
	 * Constructor.
	 *
	 * @param tau Pointer to a SingletonGraphView object with the
	 * precision parameter of the random effects as the sampled node.
	 *
	 * @param glmmethod Pointer to a GLMMethod object that updates the
	 * random effects (and other nodes that constitute a linear
	 * predictor for a GLM).
	 */
	REMethod2(SingletonGraphView const *tau,
		  GLMMethod const *glmmethod);
	~REMethod2();
	/**
	 * Updates the random effects
	 *
	 * @param rng Random number generator used for sampling
	 */
	virtual void updateSigma(RNG *rng) = 0;
	virtual void updateTau(RNG *rng) = 0;
	void update(RNG *rng);
	void calDesignSigma();
	/** 
	 * The likelihood for the standard deviation parameters sigma
	 * may be expressed in canonical form as
	 * <pre>
	 * - t(delta) %*% A %*% delta/2 + t(delta) %*% b
	 * </pre>
	 * where
	 * <pre>
	 * delta = sigma - sigma0
	 * </pre>
	 * and sigma0 is the current value.
	 *
	 * This function adds the likelihood contributions to A and b.
	 *
	 * @param A pointer to an m x m matrix. On entry, A be non-zero
	 * (representing contributions from the prior): the contributions
	 * from the likelihood are added.
	 *
	 * @param b pointer to an m-vector. On entry, b may be non-zero:
	 * the contributions from the likelihood are added.
	 *
	 * @param sigma0 current value of the standard deviation parameters.
	 * The log likelihood is standardized so that the value at sigma0
	 * is zero.
	 *
	 * @param m length of sigma vector
	 *
	 */
	void calCoefSigma(double *A, double *b, double const *sigma0,
			  unsigned int m) const;
	/**
	 * Calculates the likelihood for sigma, the vector of standard
	 * deviation parameters for the random effects
	 *
	 * @param sigma proposed value of the standard deviation parameters
	 *
	 * @param sigma0 current value of the standard deviation parameters
	 *
	 * @param m length of sigma vector
	 *
	 * @see calCoefSigma
	 */
	double logLikelihoodSigma(double const *sigma, double const *sigma0,
				  unsigned int m) const;
    };

}}

#endif /* RE_METHOD2_H_ */
