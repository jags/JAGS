#include <config.h>

#include "DensityTotalStat.h"

#include <graph/Node.h>

#include <cmath>
#include <stdexcept>

using std::vector;
using std::logic_error;
using std::exp;

namespace jags {
    namespace diag {
	
	DensityTotalStat::DensityTotalStat(vector<Node const *> const &nodes, 
					   DensityType const density_type)
	    : MonitorStat(1UL), _nodes(nodes), _density_type(density_type)
	{
	    if (density_type == DTUNSET) {
		throw logic_error("Unimplemented DensityType in DensityTotal");
	    }
	}
	
	vector<double> DensityTotalStat::value(unsigned int ch) const
	{
	    double logdensity = 0.0, loglik = 0.0;
	    for (auto p = _nodes.begin(); p != _nodes.end(); ++p) {
		double ld = (*p)->logDensity(ch, PDF_FULL);
		logdensity += ld;
		loglik += (*p)->isFixed() ? ld : 0.0;
	    }

	    double ans = 0.0;
	    switch(_density_type) {
	    case DENSITY:
		ans = exp(logdensity);
		break;
	    case LOGDENSITY:
		ans = logdensity;
		break;
	    case LIKELIHOOD:
		ans = exp(loglik);
		break;
	    case LOGLIKELIHOOD:
		ans = loglik;
		break;
	    case DEVIANCE:
		ans = -2.0 * loglik;
		break;
	    case DTUNSET:
		break;
	    }
	    return vector<double>(1, ans);
	}
	
    }

}
