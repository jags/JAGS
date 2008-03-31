#include <config.h>
#include <util/nainf.h>
#include "KLBern.h"

#include <cmath>

using std::vector;
using std::log;

#define PROB0 (*par0[0])
#define PROB1 (*par1[0])

namespace dic {

    double KLBern::divergence(vector<double const *> const &par0,
			     vector<double const *> const &par1) const
    {
	return (PROB0 * log (PROB0/PROB1) +
		(1 - PROB0) * log((1 - PROB0)/(1 - PROB1)));
    }

}
