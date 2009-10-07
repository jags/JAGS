#include <config.h>

#include "ConjugateLM.h"
#include <graph/StochasticNode.h>
#include <sampler/Updater.h>

using std::vector;
using std::string;

namespace glm {

    ConjugateLM::ConjugateLM(Updater const *updater, unsigned int chain)
	: GLMMethod(updater, chain, false)
    {
    }

    double ConjugateLM::getMean(unsigned int i) const
    {
	return _updater->stochasticChildren()[i]->parents()[0]->value(_chain)[0];
    }
    
    double 
    ConjugateLM::getPrecision(unsigned int i) const 
    {
	return _updater->stochasticChildren()[i]->parents()[1]->value(_chain)[0];
    }

    double ConjugateLM::getValue(unsigned int i) const 
    {
	return _updater->stochasticChildren()[i]->value(_chain)[0];
    }
    
    string ConjugateLM::name() const
    {
	return "ConjugateLM";
    }

    void ConjugateLM::initAuxiliary(RNG *rng) 
    {
	return; //Nothing to do
    }

    void ConjugateLM::updateAuxiliary(double *b, csn const *N, RNG *rng)
    {
	return; //Nothing to do
    }
}
