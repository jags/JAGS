#ifndef CONJUGATE_FACTORY_H_
#define CONJUGATE_FACTORY_H_

#include "ConjugateSampler.h"

#include <sampler/GibbsFactory.h>
#include <map>
#include <string>
class StochasticNode;


/**
 * @short Factory object for conjugate samplers
 */
class ConjugateFactory : public GibbsFactory
{
    std::map<const std::string, bool (*)(StochasticNode *,
					 Graph const &)> _func_table;
public:
    ConjugateFactory();
    /**
     * Provides a simple interface to the canSample member functions of
     * the various conjugate samplers.  The distribution of snode is found
     * and, if it has an associated conjugate sampler, the appropriate
     * canSample member function is dispatched
     */
    bool canSample(StochasticNode *snode, Graph const &graph) const;
    /**
     * Creates a ConjugateSampler for a stochastic node. 
     */
    Sampler *makeGibbsSampler(StochasticNode *snode, Graph const &graph, unsigned int chain) const;
};

#endif /* CONJUGATE_FACTORY_H_ */
