#ifndef GLM_METHOD_H_
#define GLM_METHOD_H_

#include <sampler/MutableSampleMethod.h>

#include <string>
#include <vector>

extern "C" {
#include <cholmod.h>
}

namespace jags {

    struct RNG;
    class GraphView;
    class SingletonGraphView;
    class REMethod2;
    
namespace glm {

    class Outcome;

    /**
     * @short Abstract class for sampling generalized linear models.
     *
     * GLMMethod provides a base class for sampling methods that work
     * on generalized linear models (GLMs).  Most of the update
     * machinery is provided by GLMMethod itself. 
     *
     * GLMMethod uses sparse matrix algebra provided by the 
     * CHOLMOD library. In the context of a hierarchical model,
     * mixed models appear identical to fixed-effect models except
     * that mixed models have a design matrix that is sparse.  The use
     * of CHOLMOD, along with auxiliary variable sampling,
     * allows us to handle both fixed and random effects in a
     * consistent way without needing to distinguish between them or
     * relying on asymptotic approximations.
     */
    class GLMMethod : public MutableSampleMethod {
    protected:
	GraphView const *_view;
	unsigned int _chain;
	std::vector<SingletonGraphView const *> _sub_views;
	std::vector<Outcome *> _outcomes;
	cholmod_sparse *_x;
	cholmod_factor *_factor; //???
	void symbolic();
	void calDesign() const;
    private:
	std::vector<bool> _fixed;
	unsigned int _length_max;
	unsigned _nz_prior;
	friend REMethod2;
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
	 */
	GLMMethod(GraphView const *view, 
		  std::vector<SingletonGraphView const *> const &sub_views,
		  std::vector<Outcome *> const &outcomes,
		  unsigned int chain);
	/**
	 * Virtual destructor
	 */
	virtual ~GLMMethod();
	/**
	 * Calculates the coefficients of the posterior distribution
	 * of the regression parameters. GLMMethod uses a canonical
	 * parametrization (b, A) such that "A" is the posterior
	 * precision of the parameters and the posterior mean "mu"
	 * solves (A %*% mu = b).
	 *
	 * @param b Dense vector such that (b = A %*% mu), where "mu"
	 * is the posterior mean and "A" is the posterior precision.
	 *
	 * @param A Posterior precision represented as a sparse matrix.
	 */
	void calCoef(double *&b, cholmod_sparse *&A);
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
	 * Returns the name of the sampler
	 */
	std::string name() const;
    };

}}

#endif /* GLM_METHOD_H_ */
