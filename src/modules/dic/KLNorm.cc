#include <config.h>
#include "KLNorm.h"
#include <cmath>

using std::vector;
using std::log;

#define mu0  (*par0[0])
#define mu1  (*par1[0])
#define tau0 (*par0[1])
#define tau1 (*par1[1])

namespace dic {

    double KLNorm::divergence(vector<double const *> const &par0,
			      vector<double const *> const &par1) const
    {
	
	return  (mu0 - mu1) * (mu0 - mu1) * (tau1 + tau0)/2 + 
	    (tau1 - tau0) * (tau1 - tau0) / (2 * tau1 * tau0);
	
    }

}
