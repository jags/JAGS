#ifndef CONJUGATE_H_
#define CONJUGATE_H_

#include <sampler/GibbsSampler.h>

class StochasticNode;
class LogicalNode;

enum ConjugateDist {
  BERN, BETA, BIN, CAT, CHISQ, DEXP, DIRCH, EXP, GAMMA,
  LNORM, LOGIS, MNORM, MULTI, NEGBIN, NORM, PAR, POIS,
  T, UNIF, WEIB, WISH, OTHERDIST
};

/*
enum ConjugateOp {
  ADD, SUBTRACT, NEG, MULTIPLY, DIVIDE, INPROD, OTHEROP
};
*/

/**
 * Convenience function that converts the name of a distribution into
 * a member of the ConjugateDist enumeration. OTHERDIST is returned
 * if the name cannot be found in the internal table.
 */
ConjugateDist getDist(StochasticNode const *snode);
/** 
 * Convenience function that converts the name of a distribution into
 * a member of the ConjugateOp enumeration. OTHEROP is returned
 * if the name cannot be found in the internal table.
 */
//ConjugateOp getOp(LogicalNode const *lnode);

class ConjugateSampler : public GibbsSampler
{
 protected:
  const ConjugateDist _target_dist;
  std::vector<ConjugateDist> _child_dist;
 public:
  ConjugateSampler(StochasticNode *snode, Graph const &graph, unsigned int chain);
  /**
   * Conjugate samplers do not have a burnin mode. This function does nothing.
   */
  void burninOff();
};

#endif /* CONJUGATE_H_ */
