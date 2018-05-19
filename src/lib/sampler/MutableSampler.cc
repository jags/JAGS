#include <sampler/MutableSampler.h>
#include <sampler/MutableSampleMethod.h>
#include <graph/StochasticNode.h>
#include <sampler/GraphView.h>

#include <stdexcept>

using std::vector;
using std::logic_error;
using std::string;

namespace jags {

    MutableSampler::MutableSampler(GraphView *gv,
				   vector<MutableSampleMethod *> const &methods,
				   std::string const &name)
	: Sampler(gv), _methods(methods), _name(name)
    {
    }

    MutableSampler::~MutableSampler()
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    delete _methods[ch];
	}
    }

    void MutableSampler::update(unsigned int ch, RNG * rng)
    {
	_methods[ch]->update(rng);
    }

    void MutableSampler::adaptOff()
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    _methods[ch]->adaptOff();
	}
    }

    bool MutableSampler::checkAdaptation() const
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    if (!_methods[ch]->checkAdaptation()) return false;
	}
	return true;
    }


    bool MutableSampler::isAdaptive() const
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    if (_methods[ch]->isAdaptive())
		return true;
	}
	return false;
    }

    string MutableSampler::name() const
    {
	return _name;
    }

} //namespace jags
