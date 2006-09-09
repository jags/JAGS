#include <config.h>
#include <distribution/DistDiscrete.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <rng/RNG.h>
#include "FiniteSampler.h"

#include <cmath>
#include <cfloat>
#include <string>
#include <stdexcept>

using std::logic_error;
using std::string;

FiniteSampler::FiniteSampler(StochasticNode *node, Graph const &graph, unsigned int chain)
  : GibbsSampler(node, graph, chain)
{
  if (!canSample(node, graph))
    throw logic_error("Attempt to construct invalid FiniteSampler");
}

FiniteSampler::~FiniteSampler()
{}

void FiniteSampler::update(RNG *rng)
{
  Distribution const *dist = node()->distribution();

  long lower = static_cast<long>(dist->lowerSupport(0,
                                                    node()->parameters(chain())));
  long upper = static_cast<long>(dist->upperSupport(0,
                                                    node()->parameters(chain())));

  long size = upper - lower + 1;
  double *lik = new double[size];
  double liksum = 0.0;
  for (long i = 0; i < size; i++) {
    double ivalue = lower + i;
    setValue(&ivalue, 1);
    lik[i] = exp(logFullConditional());
    liksum += lik[i];
  }

  /* Sample */
  double urand = rng->uniform() * liksum;
  long i;
  liksum = 0.0;
  for (i = 0; i < size - 1; i++) {
    liksum += lik[i];
    if (liksum > urand) {
      break;
    }
  }
  double ivalue = lower + i;
  setValue(&ivalue, 1);
  delete [] lik;
}

bool FiniteSampler::canSample(StochasticNode const *node,
			      Graph const &graph)
{
  //Node must be scalar with discrete-valued distribution of full rank
  Distribution const *dist = node->distribution();
  if (!dist->isDiscreteValued())
    return false;

  if (node->length() != 1)
    return false;

  if (dist->isDeterministic())
    return false;

  for (unsigned int ch = 0; ch < node->nchain(); ++ch) {

    //Distribution cannot be unbounded
    if (dist->upperSupport(0, node->parameters(ch)) == DBL_MAX ||
	dist->lowerSupport(0, node->parameters(ch)) == -DBL_MAX) {
      return false;
    }
    else {
      double n = dist->upperSupport(0, node->parameters(ch)) -
	dist->lowerSupport(0, node->parameters(ch)) + 1;
      if (n <= 1 || n > 20) //fixme: totally arbitrary
         return false;
      //fixme: should add condition that discrete-valued parents are fixed
    }
  }
  return true;
}

void FiniteSampler::burninOff()
{}
