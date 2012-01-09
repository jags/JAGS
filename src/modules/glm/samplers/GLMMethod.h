#ifndef GLM_METHOD_H_
#define GLM_METHOD_H_

#include <sampler/SampleMethod.h>
#include <sampler/GraphView.h>

#include <string>
#include <vector>

extern "C" {
#include <cholmod.h>
}

class Graph;
class RNG;
class StochasticNode;

enum GLMFamily {GLM_NORMAL, GLM_BERNOULLI, GLM_BINOMIAL, GLM_POISSON, 
		GLM_UNKNOWN};

namespace glm {

    /**
     * @short Abstract class for sampling generalized linear models.
     *
     * GLMMethod provides a base class for sampling methods that work
     * on generalized linear models (GLMs).  Most of the update
     * machinery is provided by GLMMethod itself. Sub-classes have to
     * define member functions that define the inputs to the sampling
     * algorithm.
     *
     * Classes inheriting from GLMMethod typically use auxiliary
     * variable sampling.  This is a form of data augmentation in
     * which additional latent (auxiliary) variables are added to the
     * model to reduce it to the form of a linear model with normal
     * outcome variables.  We also assume that that the regression
     * parameters have a prior multivariate normal distribution.  Thus
     * the posterior distribution is also multivariate normal and the
     * the regression parameters can be sampled efficiently in a block.
     *
     * GLMMethod uses sparse matrix algebra provided by the libraries
     * CHOLMOD and CSparse. In the context of a hierarchical model,
     * mixed models appear identical to fixed-effect models except
     * that mixed models have a design matrix that is sparse.  The use
     * of CHOLMOD and CSparse, along with auxiliary variable sampling,
     * allows us to handle both fixed and random effects in a
     * consistent way without needing to distinguish between them, or
     * rely on asymptotic approximations.
     */
    class GLMMethod : public SampleMethod {
	std::vector<double const *> _lp;
    protected:
	GraphView const *_view;
	unsigned int _chain;
	std::vector<GraphView const *> _sub_views;
	cholmod_sparse *_x;
	cholmod_factor *_factor;
	void symbolic();
	void calDesign() const;
    private:
	std::vector<bool> _fixed;
	unsigned int _length_max;
	unsigned _nz_prior;
	bool _init;
    public:
	/**
	 * Constructor.
	 *
	 * @param view Pointer to a GraphView object for all sampled nodes.
	 *
	 * @param sub_views Vector of pointers to GraphView objects with
	 * length equal to the number of sampled nodes. Each sub-view
	 * corresponds to a single sampled node.
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
	GLMMethod(GraphView const *view, 
		  std::vector<GraphView const *> const &sub_views,
		  unsigned int chain, bool link);
	/**
	 * Virtual destructor
	 */
	virtual ~GLMMethod();
	/**
	 * Updates the regression parameters by treating the GLM as a
	 * linear model (LM), either by a linear approximation or by
	 * using auxiliary variables. All regression parameters are
	 * updated together in a block.
	 *
	 * The updateLM function relies on virtual functions which are
	 * implemented by sub-classes: getPrecision, getValue and
	 * updateAuxiliary.
	 *
	 * In order to provide more flexibility, updateLM has an optional
	 * arguments stochastic.
	 *
	 * @param rng Random number generator used for sampling
	 *
	 * @param T Temperature parameter that allows a tempered (T < 1)
	 * or annealed (T > 1) update. See calCoef.
	 */
	void updateLM(RNG *rng, double T = 1);
	/**
	 * Updates the regression parameters element-wise (i.e. with
	 * Gibbs sampling).  Although Gibbs sampling less efficient
	 * than the block-updating provided by updateLM, it is more
	 * widely applicable.  In particular, if the regression
	 * parameters have a truncated normal prior, the model is
	 * still amenable to Gibbs sampling.
	 */
	void updateLMGibbs(RNG *rng);
	/**
	 * Calculates the coefficients of the posterior distribution
	 * of the regression parameters. GLMMethod uses a canonical
	 * parametrization (b, A) such that "A" is the posterior
	 * precision of the parameters and the posterior mean "mu"
	 * solves (A %*% mu = b). Both b and A are derived from the
	 * sum of prior and likelihood contributions. The calCoef
	 * function optionally allows these contributions to be
	 * weighted.
	 *
	 * @param b Dense vector such that (b = A %*% mu), where "mu"
	 * is the posterior mean and "A" is the posterior precision.
	 *
	 * @param A Posterior precision represented as a sparse matrix.
	 *
	 * @param Temp Temperature parameter that is used to weight
	 * the contribution from the likelihood.
	 */
	void calCoef(double *&b, cholmod_sparse *&A, double Temp=1.0);
	/**
	 * Returns the linear predictor for the outcome variable with index i.
	 */
	double getMean(unsigned int i) const;
	/**
	 * Returns the precision of the outcome variable with index i.
	 * This may be the precision parameter of an auxiliary variable.
	 */
	virtual double getPrecision(unsigned int i) const = 0;
	/**
	 * Returns the value of the outcome variable with index i.
	 * This may be an auxiliary variable, rather than the observed
	 * outcome.
	 */
	virtual double getValue(unsigned int i) const = 0;
	/**
	 * Updates auxiliary variables.  The default method does
	 * nothing.  Sampling methods that use auxiliary variables to
	 * reduce the GLM to a linear model must provide their own
	 * implementation.
	 *
	 * This function is called by updateLM.  Internally, updateLM
	 * calculates the posterior mean "mu" by solving the equation
	 * (A %*% mu = b) where "A" is the posterior precision. The
	 * same dense vector is used to hold "b" (before solving the
	 * equation) and "mu" (after solving the equation).  If
	 * updateLM is called with parameter "chol" set to false then
	 * updateAuxiliary is called before solving the equation: thus
	 * the first argument (y) should contain the canonical
	 * parameter (b).  If updateLM is called with "chol" set to
	 * true then updateAuxiliary is called after solving the
	 * equation: thus the first argument (y) should contain the
	 * posterior mean.
	 *
	 * IMPORTANT NOTE: GLMMethod uses a parameterization in which
	 * the current value of the parameters is considered the
	 * origin. The value of "y" (mu or b) may need to be adjusted
	 * for this centring by an implementation of updateAuxiliary.
	 * 
	 * @param y Dense vector which may be either the posterior
	 * mean "mu" or (A %*% mu), where "A" is the posterior
	 * precision.
	 *
	 * @param N Cholesky factorization of the posterior precision "A".
	 *
	 * @param rng Random number generator used for sampling.
	 */
	virtual void updateAuxiliary(cholmod_dense *y, cholmod_factor *N, 
				     RNG *rng);
	/**
	 * Returns false. Sampling methods inheriting from GLMMethod
	 * are not adaptive.
	 */
	bool isAdaptive() const;
	/**
	 * Does nothing, as GLMMethod is not adaptive.
	 */
	void adaptOff();
	/**
	 * Returns true, as GLMMethod is not adaptive
	 */
	bool checkAdaptation() const;
	/**
	 * Utility function that classifies the distribution of a
	 * stochastic node into one of the classes defined by the
	 * enumeration GLMFamily.
	 */
	static GLMFamily getFamily(StochasticNode const *snode);
    };

}
#endif /* GLM_METHOD_H_ */
