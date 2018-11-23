#include <config.h>

#include "DirichletCat.h"

#include <rng/RNG.h>
#include <graph/MixtureNode.h>
#include <graph/StochasticNode.h>
#include <sampler/GraphView.h>
#include <module/ModuleError.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

#include <JRmath.h>

using std::vector;
using std::set;
using std::string;
using std::copy;
using std::map;
using std::pair;


namespace jags {

    static vector<MixtureNode const *>  makeMixParents(GraphView const *gv) 
    {
	vector<StochasticNode *> const &schildren =
	    gv->stochasticChildren();
	
	vector<MixtureNode const *> mixnodes;
	mixnodes.reserve(schildren.size());
	for (unsigned int i = 0; i < schildren.size(); ++i) {
	    mixnodes.push_back(asMixture(schildren[i]->parents()[0]));
	}
	return mixnodes;
    }

    static map<Node const *, vector<double> > makeParMap(GraphView const *gv)
    {
	vector<StochasticNode *> const &snodes =  gv->nodes();
	map<Node const *, vector<double> > parmap;
	unsigned int N = snodes[0]->length();
	for (unsigned int i = 0; i < snodes.size(); ++i) {
	    if (snodes[i]->length() != N) {
		throwLogicError("Length mismatch in MixDirich");
	    }
	    parmap[snodes[i]] = vector<double>(N);
	}
	return parmap;
    }

    namespace mix {
	
	bool DirichletCat::canSample(GraphView const *gv)
	{
	    // Check sample nodes
	    vector<StochasticNode*> snodes = gv->nodes();
	    for (unsigned int i = 0; i < snodes.size(); ++i) {
		if (isBounded(snodes[i])) return false;
		if (snodes[i]->distribution()->name() != "ddirch")
		    return false;
		//Sampled nodes must be completely unobserved
		if (isObserved(snodes[i]))
		    return false;
	    }

	    set<Node const *> sset;
	    set<Node const *> dset;
	    for (unsigned int i = 0; i < snodes.size(); ++i) {
		sset.insert(snodes[i]);
		dset.insert(snodes[i]);
	    }

	    vector<DeterministicNode*> const &dchildren = 
		gv->deterministicChildren();
	    
	    for (unsigned int i = 0; i < dchildren.size(); ++i) {
		if (MixtureNode const *m = asMixture(dchildren[i])) {
		    vector<Node const*> const &mparents = m->parents();
		    unsigned int Nindex = m->index_size();
		    for (unsigned int j = 0; j < Nindex; ++j) {
			if (dset.count(mparents[j])) return false;
		    }
		    for (unsigned int j = Nindex; j < mparents.size(); ++j) {
			if (sset.count(mparents[j]) == 0) return false;
		    }
		}
		else {
		    return false;
		}
		dset.insert(dchildren[i]);
	    }
    
	    // Check stochastic children
	    vector<StochasticNode *> const &schildren = 
		gv->stochasticChildren();
	    
	    for (unsigned int i = 0; i < schildren.size(); ++i) {
		if (isBounded(schildren[i])) 
		    return false;
		if (schildren[i]->distribution()->name() != "dcat")
		    return false;
		if (sset.count(schildren[i]->parents()[0]))
		    return false;
	    }

	    return true;
	}
	
	DirichletCat::DirichletCat(GraphView const *gv, unsigned int chain)
	    : _gv(gv), _parmap(makeParMap(gv)), _mixparents(makeMixParents(gv)),
	      _chain(chain), _size(gv->nodes()[0]->length())
	{
	    updateParMap();
	}

	vector<double> &DirichletCat::getActiveParameter(unsigned int i)
	{
	    Node const *active = _mixparents[i]->activeParent(_chain);
	    map<Node const *, vector<double> >::iterator p =
		_parmap.find(active);
	    if (p == _parmap.end()) {
		throwLogicError("No active parameter in DirichletCat sampler");
	    }
	    return p->second;
	}
	
	void DirichletCat::updateParMap()
	{
	    vector<StochasticNode*> const &snodes = _gv->nodes();

	    //Priors
	    for (unsigned int i = 0; i < snodes.size(); ++i) {
		double const *prior = snodes[i]->parents()[0]->value(_chain);
		vector<double> &param = _parmap[snodes[i]];
		copy(prior, prior + _size, param.begin());
	    }

	    vector<StochasticNode *> const &schildren = 
		_gv->stochasticChildren();

	    //Likelihoods
	    for (unsigned int j = 0; j < schildren.size(); ++j) {
		Node const *active = _mixparents[j]->activeParent(_chain);
		vector<double> &param = _parmap.find(active)->second;
		int index = static_cast<int>(*schildren[j]->value(_chain)) - 1;
		if (param[index] > 0) {
		    param[index] += 1;
		}
		else {
		    throwNodeError(schildren[j], 
				   "Likelihood error in DirichletCat sample");
		}

	    }
	}
	
	void DirichletCat::update(RNG *rng)
	{
	    updateParMap();

	    /* 
	       Draw Dirichlet sample by drawing independent gamma random
	       variates and then normalizing
	    */

	    vector<StochasticNode*> const &snodes = _gv->nodes();
	    vector<double> xnew(_gv->length());

	    for (unsigned int i = 0; i < snodes.size(); ++i) {
		
		vector<double> const &alpha = _parmap[snodes[i]];

		double xsum = 0.0;
		unsigned int offset = i * _size;
		for (unsigned int j = 0; j < _size; ++j) {
		    if (alpha[j] > 0) {
			xnew[j + offset] = rgamma(alpha[j], 1, rng);
			xsum += xnew[j + offset];
		    }
		    else {
			xnew[j + offset] = 0;
		    }
		} 
		for (unsigned int j = 0; j < _size; ++j) {
		    xnew[j + offset] /= xsum;
		}
	    }
	    
	    _gv->setValue(xnew, _chain);
	}

    }
}

