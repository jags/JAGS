#include "DensityEnums.h"

using std::string;

namespace jags {
namespace dic {

	bool getMonitorDensityTypes(string const &type, MonitorType &monitor_type, 
		DensityType &density_type)
	{

		if (type == "density_trace") {
			monitor_type = TRACE;
			density_type = DENSITY;
		}
		else if (type == "density_mean") {
			monitor_type = MEAN;
			density_type = DENSITY;
		}
		else if (type == "density_variance") {
			monitor_type = VARIANCE;
			density_type = DENSITY;
		}
		else if (type == "density_total") {
			monitor_type = TOTAL;
			density_type = DENSITY;
		}
		else if (type == "density_poolmean") {
			monitor_type = POOLMEAN;
			density_type = DENSITY;
		}
		else if (type == "density_poolvariance") {
			monitor_type = POOLVARIANCE;
			density_type = DENSITY;
		}
		else if (type == "logdensity_trace") {
			monitor_type = TRACE;
			density_type = LOGDENSITY;
		}
		else if (type == "logdensity_mean") {
			monitor_type = MEAN;
			density_type = LOGDENSITY;
		}
		else if (type == "logdensity_variance") {
			monitor_type = VARIANCE;
			density_type = LOGDENSITY;
		}
		else if (type == "logdensity_total") {
			monitor_type = TOTAL;
			density_type = LOGDENSITY;
		}
		else if (type == "logdensity_poolmean") {
			monitor_type = POOLMEAN;
			density_type = LOGDENSITY;
		}
		else if (type == "logdensity_poolvariance") {
			monitor_type = POOLVARIANCE;
			density_type = LOGDENSITY;
		}
		else if (type == "deviance_trace") {
			monitor_type = TRACE;
			density_type = DEVIANCE;
		}
		else if (type == "deviance_mean") {
			monitor_type = MEAN;
			density_type = DEVIANCE;
		}
		else if (type == "deviance_variance") {
			monitor_type = VARIANCE;
			density_type = DEVIANCE;
		}
		else if (type == "deviance_total") {
			monitor_type = TOTAL;
			density_type = DEVIANCE;
		}
		else if (type == "deviance_poolmean") {
			monitor_type = POOLMEAN;
			density_type = DEVIANCE;
		}
		else if (type == "deviance_poolvariance") {
			monitor_type = POOLVARIANCE;
			density_type = DEVIANCE;
		}
		else if (type == "pD") {
			monitor_type = PD;
		}
		else if (type == "pD_total") {
			monitor_type = PDTOTAL;
		}
		else if (type == "popt") {
			monitor_type = POPT;
		}
		else if (type == "popt_total") {
			monitor_type = POPTTOTAL;
		}
		else if (type == "popt_total_replicatemean") {
			monitor_type = POPTTOTALREP;
		}
		else if (type == "pv") {
			monitor_type = PV;
		}
		else {
			// If not listed above:
			return false;
		}
	
		return true;
		
	}
 
}}
