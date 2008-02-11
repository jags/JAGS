#ifndef CONJUGATE_SAMPLER_H_
#define CONJUGATE_SAMPLER_H_

#include <sampler/Sampler.h>
#include "ConjugateMethod.h"

class StochasticNode;
class LogicalNode;

enum ConjugateDist {
    BERN, BETA, BIN, CAT, CHISQ, DEXP, DIRCH, EXP, GAMMA,
    LNORM, LOGIS, MNORM, MULTI, NEGBIN, NORM, PAR, POIS,
    T, UNIF, WEIB, WISH, OTHERDIST
};

/**
 * Convenience function that converts the name of a distribution into
 * a member of the ConjugateDist enumeration. OTHERDIST is returned
 * if the name cannot be found in the internal table.
 */
ConjugateDist getDist(StochasticNode const *snode);

class ConjugateSampler : public Sampler
{
    StochasticNode *_snode;
    ConjugateMethod *_method;
    const ConjugateDist _target_dist;
    std::vector<ConjugateDist> _child_dist;
public:
    ConjugateSampler(StochasticNode* node, Graph const &graph, 
		     ConjugateMethod *method);
    ~ConjugateSampler();
    StochasticNode *node() const;
    ConjugateDist targetDist() const;
    std::vector<ConjugateDist> const &childDist() const;
    void update(std::vector<RNG*> const &);
    /**
     * Conjugate methods are not adaptive
     */
    bool isAdaptive() const;
    bool adaptOff();
    std::string name() const;
};

#endif /* CONJUGATE_SAMPLER_H_ */
