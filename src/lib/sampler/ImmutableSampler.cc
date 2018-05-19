#include <config.h>
#include <sampler/ImmutableSampler.h>
//Needed for name
#include <sampler/ImmutableSampleMethod.h>
//Needed for nchain
#include <sampler/GraphView.h>

using std::vector;
using std::string;

namespace jags {

    ImmutableSampler::ImmutableSampler(GraphView *gv, 
				       ImmutableSampleMethod *method,
				       std::string const &name)
	: Sampler(gv), _method(method), _nchain(nchain(gv)), _name(name)
    {
    }

    ImmutableSampler::~ImmutableSampler()
    {
	delete _method;
    }

    void ImmutableSampler::update(unsigned int ch, RNG * rng)
    {
	_method->update(ch, rng);
    }

    bool ImmutableSampler::isAdaptive() const
    {
	return false;
    }

    void ImmutableSampler::adaptOff()
    {
    }

    bool ImmutableSampler::checkAdaptation() const
    {
	return true;
    }

    string ImmutableSampler::name() const
    {
	return _name;
    }

}
