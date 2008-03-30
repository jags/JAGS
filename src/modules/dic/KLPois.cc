#include <config.h>
#include "KLPois.h"

#include <cmath>

using std::vector;
using std::log;

#define lambda0  (*par0[0])
#define lambda1  (*par1[0])

namespace dic {

    double KLPois::divergence(vector<double const *> const &par0,
			      vector<double const *> const &par1) const
    {
	
	return  lambda0 * log(lambda0/lambda1) + (lambda1 - lambda0);
    }
    
}
