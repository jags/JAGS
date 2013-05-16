#ifndef GLM_OUTCOME_H_
#define GLM_OUTCOME_H_

namespace jags {

class RNG;
class StochasticNode;

namespace glm {
    
    /**
     * @short Outcome for a generalized linear model.
     *
     * The Outcome class provides an abstract interface to the outcome
     * variable of a generalized linear model. Some GLMs can be
     * reduced to a normal linear model by augmentation with auxiliary
     * (latent) variables. Other GLM methods use a local linear
     * approximation.  The Outcome class encapsulates this by
     * providing member functions "value", "mean" and "precision",
     * corresponding to the normal approximation, as well as functions
     * for updating any auxiliary variables, which may be hidden
     * inside the Outcome object.
     */
    class Outcome  {
      protected:
	double const &_lp;
      public:
	/**
	 * Constructor
	 *
	 * @param snode Stochastic node representing the true outcome
	 * variable of a GLM 
	 * 
	 * @param chain Index number of the chain (starting from zero)
	 * to use
	 */
	Outcome(StochasticNode const *snode, unsigned int chain);
	virtual ~Outcome();
	/**
	 * Returns the value of the linear predictor, which is the
	 * mean of the normal approximation.
	 */
	virtual double mean() const;
	/**
	 * Returns the precision of the normal approximation.
	 */
	virtual double precision() const = 0;
	/**
	 * Returns the value of the normal approximation.
	 */
	virtual double value() const = 0;
	/**
	 * Updates the auxiliary variables using the current value of
	 * the linear predictor. The default implementation does
	 * nothing.
	 *
	 * @param rng Random number generator
	 */
	virtual void update(RNG *rng); 
	/**
	 * Updates the auxiliary variables marginalizing over the
	 * distribution of the linear predictor. The default
	 * implementation does nothing.
	 *
	 * @param mean Expected value of the linear predictor.
	 * @param var  Variance of the linear predictor.
	 * @param rng  Random number generator.
	 */
	//virtual void update(double mean, double var, RNG *rng);
    };

}}

#endif /* GLM_OUTCOME_H_ */