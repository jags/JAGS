#include <config.h>

#include "DensityStat.h"

#include <graph/Node.h>

#include <cmath>
#include <stdexcept>

using std::vector;
using std::exp;
using std::logic_error;

#include <graph/Node.h>

namespace jags {
    namespace diag {


	static inline double ld(Node const *node, unsigned int ch) {
	    return node->logDensity(ch, PDF_FULL);
	}
	
	static inline double ll(Node const *node, unsigned int ch) {
	    return node->logLikelihood(ch);
	}
	
	DensityStat::DensityStat(vector<Node const *> const &nodes, 
				 DensityType density_type)
	    : MonitorStat(nodes.size()), _nodes(nodes), _density_type(density_type)
	{
	    if (density_type == DTUNSET) {
		throw logic_error("Unimplemented DensityType in DensityStat");
	    }
	}

	vector<double> DensityStat::value(unsigned int ch) const
	{
	    vector<double> value(_nodes.size());
	    for (unsigned long i = 0; i < _nodes.size(); ++i) {
		switch(_density_type) {
		case DENSITY:
		    value[i] = exp(ld(_nodes[i], ch));
		    break;
		case LOGDENSITY:
		    value[i] = ld(_nodes[i], ch);
		    break;
		case LIKELIHOOD:
		    value[i] = exp(ll(_nodes[i], ch));
		    break;
		case LOGLIKELIHOOD:
		    value[i] = ll(_nodes[i], ch);
		    break;
		case DEVIANCE:
		    value[i] = -2.0 * ll(_nodes[i], ch);
		    break;
		case DTUNSET:
		    break; //-Wswitch
		}
	    }
	    return value;
	}
		
    }
}
