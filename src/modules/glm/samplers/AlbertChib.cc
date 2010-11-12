#include <config.h>

#include "AlbertChib.h"
#include "KS.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <sampler/GraphView.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <cmath>

#include <stdexcept>

using std::vector;
using std::exp;
using std::fabs;
using std::string;
using std::logic_error;
using std::log;

//FIXME: maybe use R math library here

//Left truncated logit
static double llogit(double left, RNG *rng, double mu)
{
    double qleft = 1/(1 + exp(mu-left));
    double x = qleft + (1 - qleft) * rng->uniform();
    return mu + log(x) - log(1 - x);
}

//Right truncated logit
static double rlogit(double right, RNG *rng, double mu)
{
    double qright = 1/(1 + exp(mu-right));
    double x = qright * rng->uniform();
    return mu + log(x) - log(1 - x);
}

#define CHILD(i) (_view->stochasticChildren()[i])

namespace glm {

    AlbertChib::AlbertChib(GraphView const *view,
			   vector<GraphView const *> const &sub_views,
			   unsigned int chain, bool gibbs)
	: BinaryGLM(view, sub_views, chain), _gibbs(gibbs)
    {
    }

    void AlbertChib::update(RNG *rng)
    {

	if (_gibbs) {
	    updateLMGibbs(rng);
	}
	else {
	    updateLM(rng);
	}

	unsigned int nrow = _view->stochasticChildren().size();
	
	double y, mu;
	for (unsigned int r = 0; r < nrow; ++r) {
	    switch(_outcome[r]) {
	    case BGLM_NORMAL:
		break;
	    case BGLM_PROBIT:
		y = CHILD(r)->value(_chain)[0];
		if (y == 1) {
		    _z[r] = lnormal(0, rng, getMean(r));
		}
		else if (y == 0) {
		    _z[r] = rnormal(0, rng, getMean(r));
		}
		else {
		    throw logic_error("Invalid child value in HolmesHeld");
		}
		break;
	    case BGLM_LOGIT:
		y = CHILD(r)->value(_chain)[0];
		mu = getMean(r);
		if (y == 1) {
		    _z[r] = llogit(0, rng, mu);
		}
		else if (y == 0) {
		    _z[r] = rlogit(0, rng, mu);
		}
		else {
		    throw logic_error("Invalid child value in HolmesHeld");
		}
		_tau[r] = 1/sample_lambda(fabs(_z[r] - mu), rng);
		break;
	    }
	}


	
    }
	
    string AlbertChib::name() const
    {
	if (_gibbs)
	    return "Albert-Chib-Gibbs";
	else
	    return "Albert-Chib";
    }
}
    
