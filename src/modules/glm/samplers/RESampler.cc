#include "RESampler.h"
#include "REMethod.h"

#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/GraphView.h>

using std::vector;
using std::string;

namespace jags {
    namespace glm {
	
	RESampler::RESampler(GraphView *view,
			     SingletonGraphView *tau, GraphView *eps,
			     vector<SingletonGraphView*> sub_eps,
			     vector<REMethod*> const &methods,
			     string const &name)
	    : Sampler(view), _tau(tau), _eps(eps), _sub_eps(sub_eps),
	      _methods(methods), _name(name)
	{
	}

	RESampler::~RESampler()
	{
	    delete _tau;
	    delete _eps;
	    for (unsigned int i = 0; i < _sub_eps.size(); ++i) {
		delete _sub_eps[i];
	    }
	    for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
		delete _methods[ch];
	    }
	}

	void RESampler::update(vector<RNG*> const &rngs)
	{
	    for (unsigned int ch = 0; ch < rngs.size(); ++ch) {
		_methods[ch]->update(rngs[ch]);
	    }
	}

	void RESampler::adaptOff()
	{
	    for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
		_methods[ch]->adaptOff();
	    }
	}

	bool RESampler::checkAdaptation() const
	{
	    for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
		if (!_methods[ch]->checkAdaptation()) return false;
	    }
	    return true;
	}


	bool RESampler::isAdaptive() const
	{
	    for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
		if (_methods[ch]->isAdaptive())
		    return true;
	    }
	    return false;
	}

	string RESampler::name() const
	{
	    return _name;
	}

    } //namespace glm
} //namespace jags
