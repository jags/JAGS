#include <config.h>
#include <util/nainf.h>
#include "KLBin.h"

#include <cmath>

using std::vector;
using std::log;

#define PROB0 (*par0[0])
#define PROB1 (*par1[0])
#define N0 (*par0[1])
#define N1 (*par1[1])

namespace dic {

    double KLBin::divergence(vector<double const *> const &par0,
			     vector<double const *> const &par1) const
    {
	if (N0 != N1) {
	    return JAGS_POSINF;
	}
	else {
	    return N0 * ((PROB0 - PROB1) * 
			 (log (PROB0/(1 - PROB0)) - log(PROB1/(1 - PROB1))));
	}
    }

}
