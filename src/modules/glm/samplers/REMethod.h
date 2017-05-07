#ifndef RE_METHOD_H_
#define RE_METHOD_H_

#include "GLMMethod.h"

#include <vector>

namespace jags {

    struct RNG;
    class GraphView;
    class SingletonGraphView;

namespace glm {

    class Outcome;

    /**
     * @short Sample method for random effects and their variance.
     */
    class REMethod : public GLMMethod {
      protected:
	SingletonGraphView const *_tau;
	GraphView const *_eps;
	cholmod_dense *_z;
      public:
	/**
	 * Constructor.
	 *
	 * @param view Pointer to a GraphView object for all sampled nodes.
	 *
	 * @param sub_views Vector of pointers to SingletonGraphView
	 * objects with length equal to the number of sampled
	 * nodes. Each sub-view corresponds to a single sampled node.
	 *
	 * @param outcomes Vector of pointers to Outcome objects with length
	 * equal to the number of stochastic children of the sampled nodes.
	 * The GLMMethod objects takes ownership of each Outcome in the vector
	 * and frees the memory in the destructor.
	 * 
	 * @param chain Number of the chain (starting from 0) to which
	 * the sampling method will be applied.
	 *
	 * @param link Boolean flag that is passed to the utility
	 * function checkLinear when checking to see if we have a
	 * valid GLM. If link is true then the last deterministic
	 * descendents in view (i.e. those with no deterministic
	 * descendants) may be link nodes.
	 */
	REMethod(SingletonGraphView const *tau,
		 GraphView const *eps, 
		 std::vector<SingletonGraphView const *> const &sub_eps,
		 std::vector<Outcome *> const &outcomes,
		 unsigned int chain);
	~REMethod();
	/**
	 * Updates the random effects
	 *
	 * @param rng Random number generator used for sampling
	 */
	void updateEps(RNG *rng);
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
	/**
	 * Rescales random effects after sampling sigma
	 *
	 * @param sigma proposed value of the standard deviation parameters
	 *
	 * @param sigma0 current value of the standard deviation parameters
	 *
	 * @param m length of sigma vector
	 */
	void rescaleSigma(double const *sigma, double const *sigma0,
			  unsigned int m);
    };

}}

#endif /* RE_METHOD_H_ */
