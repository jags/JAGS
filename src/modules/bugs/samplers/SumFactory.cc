#include <config.h>
//#include <distribution/Distribution.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
//#include <graph/NodeError.h>
#include <sampler/MutableSampler.h>
#include <sampler/GraphView.h>

#include "SumFactory.h"
#include "SumMethod.h"

//#include <algorithm>
#include <map>

using std::list;
using std::vector;
using std::string;
using std::map;

namespace jags {
    namespace bugs {

	vector<Sampler *>
	SumFactory::makeSamplers(list<StochasticNode*> const &nodes,
				 Graph const &graph) const
	{
	    vector<Sampler*> samplers;
	    
	    vector<StochasticNode const*> sum_nodes;
	    map<StochasticNode const*, vector<StochasticNode*> > smap;
	    map<StochasticNode const*, vector<StochasticNode*> >::iterator i;

	    for (list<StochasticNode*>::const_iterator p = nodes.begin();
		 p != nodes.end(); ++p)
	    {
		StochasticNode const *sum = SumMethod::isCandidate(*p, graph);
		if (sum) {
		    i = smap.find(sum);
		    if (i == smap.end()) {
			sum_nodes.push_back(sum);
			smap[sum] = vector<StochasticNode*>(1, *p);
		    }
		    else {
			i->second.push_back(*p);
		    }
		}
	    }

	    for (vector<StochasticNode const *>::const_iterator p =
		     sum_nodes.begin(); p != sum_nodes.end(); ++p)
	    {
		i = smap.find(*p);
		if (!SumMethod::canSample(i->second, graph)) continue;

		GraphView *gv = new GraphView(i->second, graph, true);		
		unsigned int N = nchain(gv);
		vector<MutableSampleMethod*> methods(N);
		for (unsigned int ch = 0; ch < N; ++ch) {
		    methods[ch] = new SumMethod(gv, ch);
		}
		Sampler *s = new MutableSampler(gv, methods, "bugs::Sum");
		samplers.push_back(s);
	    }

	    return samplers;
	}

	string SumFactory::name() const
	{
	    return "bugs::Sum";
	}

    }
}
