#ifndef CLASSIFY_H_
#define CLASSIFY_H_

namespace jags {

class StochasticNode;

namespace glm {

    /**
     * Enumerates distribution families for generalized linear models
     */
    enum GLMFamily {GLM_NORMAL, GLM_BERNOULLI, GLM_BINOMIAL, GLM_POISSON,
		    GLM_LOGISTIC, GLM_T, GLM_ORDLOGIT, GLM_ORDPROBIT,
		    GLM_UNKNOWN};
    
    /**
     * Enumerates link functions for generalized linear models
     */
    enum GLMLink {LNK_LINEAR, LNK_LOG, LNK_LOGIT, LNK_PROBIT, LNK_UNKNOWN};

    /**
     * Utility function that classifies the distribution of a
     * stochastic node into one of the classes defined by the
     * enumeration GLMFamily.
     */
    GLMFamily getFamily(StochasticNode const *snode);

    /**
     * Utility function that classifies the link function of a
     * stochastic node into one of the classes defined by the
     * enumeration GLMLink.
     */
    GLMLink getLink(StochasticNode const *snode);

}}

#endif /* CLASSIFY_H_ */
