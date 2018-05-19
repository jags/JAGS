#include <config.h>

#include "GLMSampler.h"
#include "GLMMethod.h"

#include <sampler/GraphView.h>
#include <sampler/SingletonGraphView.h>


using std::vector;
using std::string;

namespace jags {
namespace glm {

    GLMSampler::GLMSampler(GraphView *view, 
			   vector<SingletonGraphView*> const &sub_views,
			   vector<GLMMethod*> const &methods,
			   std::string const &name)
	: Sampler(view), _view(view), _sub_views(sub_views), _methods(methods),
	  _name(name)
    {
	//FIXME We need a pointer _view here because the friend class
	//REFactory2 needs access to it and cannot get it from the parent
	//Sampler class. Need to add an extractor function.
    }
    
    GLMSampler::~GLMSampler()
    {
	while (!_sub_views.empty()) {
	    delete _sub_views.back();
	    _sub_views.pop_back();
	}

	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    delete _methods[ch];
	}
    }

    void GLMSampler::update(unsigned int ch, RNG * rng)
    {
	_methods[ch]->update(rng);
    }

    void GLMSampler::adaptOff()
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    _methods[ch]->adaptOff();
	}
    }

    bool GLMSampler::checkAdaptation() const
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    if (!_methods[ch]->checkAdaptation()) return false;
	}
	return true;
    }


    bool GLMSampler::isAdaptive() const
    {
	for (unsigned int ch = 0; ch < _methods.size(); ++ch) {
	    if (_methods[ch]->isAdaptive())
		return true;
	}
	return false;
    }

    string GLMSampler::name() const
    {
	return _name;
    }

    
    vector<GLMMethod*> const &GLMSampler::methods() {
	return _methods;
    }
    

}}
