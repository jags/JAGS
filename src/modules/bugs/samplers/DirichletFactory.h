#ifndef DIRICHLET_FACTORY_H_
#define DIRICHLET_FACTORY_H_

#include "Dirichlet.h"

#include <sampler/SingletonFactory.h>
class StochasticNode;

/**
 * @short Factory object for multivariate normal samplers
 */
class DirichletFactory : public SingletonFactory
{
public:
    bool canSample(StochasticNode *snode, Graph const &graph) const;
    Sampler *makeSingletonSampler(StochasticNode *snode, Graph const &graph)
                                  const;
};

#endif /* DIRICHLET_FACTORY_H_ */
