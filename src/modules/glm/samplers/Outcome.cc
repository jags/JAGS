#include "Outcome.h"
#include "Classify.h"

#include <module/ModuleError.h>
#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>

namespace jags {
namespace glm {

    static Node const *getLinearPredictor(StochasticNode const *snode)
    {
	Node const *lp = 0;

	//FIXME: Apart from checking that getFamily(snode) != GLM_UNKNOWN
	//this switch statement seems pointless.
	switch(getFamily(snode)) {
	case GLM_NORMAL: case GLM_BERNOULLI: case GLM_BINOMIAL: case GLM_POISSON:
	case GLM_LOGISTIC: case GLM_T: case GLM_ORDLOGIT:
	    lp = snode->parents()[0];
	    break;
	case GLM_UNKNOWN:
	    throwLogicError("Invalid distribution in GLMOutcome");
	    break;
	}
	
	LinkNode const *ln = dynamic_cast<LinkNode const*>(lp);
	if (ln) {
	    lp = ln->parents()[0];
	}
	
	return lp;
    }

    Outcome::Outcome(StochasticNode const *snode, unsigned int chain)
	: _lp(getLinearPredictor(snode)->value(chain)[0])
    {
    }

    Outcome::~Outcome()
    {
    }

    double Outcome::mean() const
    {
	return _lp;
    }

    void Outcome::update(RNG *rng)
    {
    }

    void Outcome::update(double mean, double var, RNG *rng)
    {
    }

    double Outcome::logMHRatio() const
    {
	return 0;
    }

    bool Outcome::fixedb() const
    {
	return false;
    }

    bool Outcome::fixedA() const
    {
	return false;
    }
}}
    

