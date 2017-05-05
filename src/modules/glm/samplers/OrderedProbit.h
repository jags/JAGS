#ifndef ORDERED_PROBIT_H_
#define ORDERED_PROBIT_H_

#include <config.h>

#include "Outcome.h"
#include <graph/StochasticNode.h>

namespace jags {

class StochasticNode;

namespace glm {

    /*
     * @short Ordered categorical outcome with probit link
     *
     * Outcome for ordered logistic regression models, based on Holmes
     * C and Held L (2006).  Bayesian Auxiliary Variables Models for
     * Binary and Multinomial Regression, Bayesian Analysis,
     * 1:148-168.
     *
     * The categorical outcome is represented in terms of a latent normal
     * normals
     */
    class OrderedProbit : public Outcome
    {
	double const &_y;
	double const *_cuts;
	const unsigned int _ncut;
	double _z;
	double _tau;
	double _sigma2;
      public:
	OrderedProbit(StochasticNode const *snode, unsigned int chain);
	double value() const;
	double precision() const;
	void update(RNG *rng);
	void update(double mean, double var, RNG *rng);
	static bool canRepresent(StochasticNode const *snode);
    };

}}

#endif /* ORDERED_PROBIT_H_ */
