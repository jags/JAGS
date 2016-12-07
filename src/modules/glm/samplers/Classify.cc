#include "Classify.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>

#include <string>

using std::vector;
using std::string;

namespace jags {
namespace glm {

    GLMFamily getFamily(StochasticNode const *snode)
    {
	string const &name = snode->distribution()->name();
	if (name == "dbern") {
	    return GLM_BERNOULLI;
	}
	else if (name == "dbin") {
	    return GLM_BINOMIAL;
	}
	else if (name == "dpois") {
	    return GLM_POISSON;
	}
	else if (name == "dnorm") {
	    return GLM_NORMAL;
	}
	else if (name == "dlogis") {
	    return GLM_LOGISTIC;
	}
	else if (name == "dt") {
	    return GLM_T;
	}
	else {
	    return GLM_UNKNOWN;
	}
    }

    GLMLink getLink(StochasticNode const *snode)
    {
	vector<Node const *> const &param = snode->parents();
	LinkNode const *lnode = dynamic_cast<LinkNode const*>(param[0]);

	if (lnode) {
	    string const &linkname = lnode->linkName();
	    if (linkname == "log")
		return LNK_LOG;
	    else if (linkname == "logit")
		return LNK_LOGIT;
	    else if (linkname == "probit")
		return LNK_PROBIT;
	    else 
		return LNK_UNKNOWN;
	}
	else {
	    return LNK_LINEAR;
	}
    }

}}
