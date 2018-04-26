#include <config.h>
#include "DirichletInfo.h"

#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <cmath>

using std::log;

namespace jags {
namespace mix {

    DirichletInfo::DirichletInfo(StochasticNode const *snode, unsigned long s,
				 unsigned int chain)
	: start(s), end(s + snode->length()), length(snode->length()),
	  sum(1), shape(0)
    {
	for (unsigned int i = 0; i < snode->length(); ++i) {
	    shape += snode->parents()[0]->value(chain)[i];
	}
    }

    double DirichletInfo::gammaPenalty() const
    {
	return (shape - 1) * log(sum) - shape * sum;
    }

}}
