#include "Outcome.h"
#include "Classify.h"

#include <module/ModuleError.h>
#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>

namespace jags {
namespace glm {

    static Node const *getLinearPredictor(StochasticNode const *snode)
    {
	if (getFamily(snode) == GLM_UNKNOWN) {
	    throwLogicError("Invalid distribution in glm::Outcome");
	}
	Node const *lp = snode->parents()[0];
	
	LinkNode const *ln = dynamic_cast<LinkNode const*>(lp);
	if (ln) {
	    lp = ln->parents()[0];
	}
	
	return lp;
    }

    Outcome::Outcome(StochasticNode const *snode, unsigned int chain)
	: _lp(getLinearPredictor(snode)->value(chain)[0]),
	  _length(snode->length()),
	  _vmean(getLinearPredictor(snode)->value(chain))
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

    unsigned int Outcome::length() const
    {
	return _length;
    }

    
    double const *Outcome::vmean() const
    {
	return _vmean;
    }
    
    double const *Outcome::vprecision() const
    {
	//FIXME: Bit of a hack 
	return 0;
    }

    double const *Outcome::vvalue() const
    {
	//FIXME: Bit of a hack 
	return 0;
    }

}}
    

