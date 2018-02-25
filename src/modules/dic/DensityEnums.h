#ifndef DENSITY_ENUMS_H_
#define DENSITY_ENUMS_H_


/**
* @short Type of density/log density/deviance calculation
*
* This enum is used by some monitors/factories in the DIC module
* to generalise calculation of values related to the deviance
*/
enum DensityType {DENSITY, LOGDENSITY, DEVIANCE};

/**
* @short Type of summarisation for deviance monitors
*
* This enum is used by some monitors/factories in the DIC module
* to generalise calculation of values related to the deviance
*/
enum MonitorType {TRACE, MEAN, VARIANCE, TOTAL, POOLMEAN, POOLVARIANCE};


#endif /* DENSITY_ENUMS_H_ */
