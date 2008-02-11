#include <config.h>
#include <Module.h>
#include <compiler/Compiler.h>
#include <function/InverseLinkFunc.h>
#include <function/Function.h>
#include <model/Model.h>
#include <algorithm>

using std::vector;
using std::list;

Module::Module()
{}

Module::~Module()
{
    unsigned int i;
    for (i = 0; i < _internal_functions.size(); ++i) {
	Compiler::funcTab().erase(_internal_functions[i]);
    }
    for (i = 0; i < _internal_link_functions.size(); ++i) {
	Compiler::funcTab().erase(_internal_link_functions[i]);
    }
    for (i = 0; i < _distributions.size(); ++i) {
	Compiler::distTab().erase(_distributions[i]);
    }

    list<RNGFactory *> &rngf = Model::rngFactories();
    for (i = 0; i < _rng_factories.size(); ++i) {
	RNGFactory *f = _rng_factories[i];
	for (;;) {
	    /* A little paranoiac, in case someone inserted the
	       same sampler factory more than once */
	    list<RNGFactory *>::iterator p = find(rngf.begin(), rngf.end(), f);
	    if (p == rngf.end())
		break;
	    rngf.erase(p);
	}
    }

    list<SamplerFactory const*> &sf = Model::samplerFactories();
    for (i = 0; i < _sampler_factories.size(); ++i) {
	SamplerFactory *f = _sampler_factories[i];
	for (;;) {
	    /* A little paranoiac, in case someone inserted the
	       same sampler factory more than once */
	    list<SamplerFactory const*>::iterator p
		= find(sf.begin(), sf.end(), f);
	    if (p == sf.end())
		break;
	    sf.erase(p);
	}
    }
}

void Module::insert(InverseLinkFunc *lfunc)
{
    _internal_link_functions.push_back(lfunc);
    _functions.push_back(lfunc);
    Compiler::funcTab().insert(lfunc);
}

void Module::insert(Function *func)
{
    _internal_functions.push_back(func);
    _functions.push_back(func);
    Compiler::funcTab().insert(func);
}

void Module::insert(Distribution *dist)
{
    _distributions.push_back(dist);
    Compiler::distTab().insert(dist);
}

void Module::insert(SamplerFactory *fac)
{
    _sampler_factories.push_back(fac);
    Model::samplerFactories().push_front(fac);
}

void Module::insert(RNGFactory *fac)
{
    _rng_factories.push_back(fac);
    Model::rngFactories().push_front(fac);
}

void Module::insert(MonitorFactory *fac)
{
    _monitor_factories.push_back(fac);
    Model::monitorFactories().push_front(fac);
}

vector<Function*> const  &Module::functions() const
{
    return _functions;
}

vector<Distribution*> const  &Module::distributions() const
{
    return _distributions;
}

vector<SamplerFactory*> const  &Module::samplerFactories() const
{
    return _sampler_factories;
}

vector<RNGFactory*> const &Module::rngFactories() const
{
    return _rng_factories;
}

vector<MonitorFactory*> const &Module::monitorFactories() const
{
    return _monitor_factories;
}
