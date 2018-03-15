#ifndef DENSITY_ENUMS_H_
#define DENSITY_ENUMS_H_

#include <string>

namespace jags {
namespace dic {

	/**
	* @short Type of density/log density/deviance calculation
	*
	* This enum is used by some monitors/factories in the DIC module
	* to generalise calculation of values related to the deviance
	*/
	enum DensityType {DTUNSET, DENSITY, LOGDENSITY, DEVIANCE};

	/**
	* @short Type of summarisation for deviance monitors
	*
	* This enum is used by some monitors/factories in the DIC module
	* to generalise calculation of values related to the deviance
	*/
	enum MonitorType {MTUNSET, TRACE, MEAN, VARIANCE, TOTAL, POOLMEAN, 
		POOLVARIANCE, PD, PDTOTAL, POPT, POPTTOTAL, POPTTOTALREP, PV};
	
	/**
	* @short Process a string to return density and monitor types
	*
	* Used by NodeDensityMonitorFactory and ObsStochDensMonitorFactory
	*/
	bool getMonitorDensityTypes(std::string const &type, 
		MonitorType &monitor_type, DensityType &density_type);
	
}
}

#endif /* DENSITY_ENUMS_H_ */
