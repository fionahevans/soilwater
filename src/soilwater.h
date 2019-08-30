#ifndef SOILWATER_H
#define SOILWATER_H

/*
	To be called by R function
*/

void soil_water_R(		
	int 	*ndays, // number of days						
	double	*rainfall, // weather data
	double	*min_temp, 									
	double	*max_temp, 									
	double	*solar_radiation,						
	int     *nlayers,  // soil parameters			
	double  *upperlimitstage1evap,   	
	double	*timeinstage2evap,			 		
	double	*evapstage1, 						 		
	double	*evapstage2, 	
	double	*diffusConst, 						 		
	double	*diffusSlope,	
	double  *runoffCurveNumber,
	double  *baresoilalbedo,					 		
	double	*time, 
	double	*tsw,
	double	*tll,
	double	*tsat,
	double  *layers_depth,
	double  *layers_airdry,	
	double  *layers_bd,
	double  *layers_wf,
	double  *layers_dryLowerLimit,
	double  *layers_drainedUpperLimit,
	double  *layers_saturation,							 		
	double  *layers_soilWaterConstant,
	double  *layers_soilWaterContent,
	double  *layers_extrSoilWaterCM,
	int     *nstages, // number of crop stages
	int     *crop_stages_in_days, // index of the start date for crop ET stages,
	double  *crop_stages_coef,  // crop coeficient parameters
	double  *soilWaterContent); // return

	
/* Calculate potential evaporation. */
double potential_evaporation(
	double min_temp, // weather data
	double max_temp,
	double solar_radiation,
	double baresoilalbedo); // soil parameter

/* 
Calculate actual evaporation for a given soil type using the potential
soil evaporation, upper limit of evaporation and amount of water infiltration.
Uses APSIM method hacked from APSIM Fortran code.
*/
double soilwat2_ritchie_evaporation(
	double  upperlimitstage1evap, // soil parameters
	double	timeinstage2evap,			 
	double	*evapstage1, 						 
	double	*evapstage2, 						
	double	*time, 									 
	double	potentialEvap, 
	double	upperLimitEvap, 
	double	waterInfiltration); // rainfall - runoff

/* 
 Calculate actual evaporation using method in HowLeakyEngine.cs.
*/

double soilwat_howleaky_evaporation(
    double  upperlimitstage1evap,   // soil parameters
    double	timeinstage2evap,			 
    double	*evapStage1, 						 
    double	*evapStage2, 						
    double	*time, 									 
    double	potentialEvap, 
    int     nlayers,
    double  *layers_soilWaterContent,
    double  *layers_airdry, 
    double	waterInfiltration); 
	
/*
Calculate runoff and remove from waterInfiltration.
*/	
double runoff(
	int    nlayers,	// soil parameters
	double runoffCurveNumber,
	double tsw,
	double tsat,
	double *layers_wf,
	double *layers_soilWaterContent,
	double *layers_dryLowerLimit,
	double *layers_extrSoilWaterCM,
	double  waterInfiltration,  // weather data
	double  previousWaterInfiltration);

/*
Calculate downward flow and redistribute water through the profile.
*/
int flux(
	int    nlayers,	// soil parameters
	double *layers_saturation,
	double *layers_soilWaterContent,
	double *layers_soilWaterConstant,
	double *layers_depth,
	double *layers_drainedUpperLimit,
	double  waterInfiltration);  // waterInfiltration
	
/*
Redistribute water between layers.  
*/	
double flow(
	int idrsw, // from flux
	int nlayers, // soil parameters
	double diffusConst, 
	double diffusSlope,
	double *layers_soilWaterContent, 
	double *layers_dryLowerLimit, 
	double *layers_depth, 
	double *layers_saturation, 
	double *layers_drainedUpperLimit,
	double tll, 
	double *tsw); 	
	
	
/* 
Iterate soil water model for one day. 

The water balance calculation method is derived from Ritchie, J. T. 
(1972). Model for predicting evaporation from a crop with incomplete 
cover. Water Resour. Res., 8, 1204?1213

FAO calculation of crop evapotranspiration from
'http://www.fao.org/docrep/x0490e/x0490e0b.htm#chapter 6'
*/
double soil_water_day_evapotranspiration(											
	double	waterInfiltration,  // weather data
	double	previousWaterInfiltration,  
	double	min_temp, 									
	double	max_temp, 									
	double	solar_radiation,						
	int     nlayers,  // soil parameters			
	double  upperlimitstage1evap,   	
	double	timeinstage2evap,			 		
	double	*evapStage1, 						 		
	double	*evapStage2, 	
	double	diffusConst, 						 		
	double	diffusSlope,	
	double  runoffCurveNumber,
	double  baresoilalbedo,					 		
	double	*time, 
	double	*tsw,
	double	tll,
	double	tsat,
	double  *layers_depth,
	double  *layers_airdry,	
	double  *layers_bd,
	double  *layers_wf,
	double  *layers_dryLowerLimit,
	double  *layers_drainedUpperLimit,
	double  *layers_saturation,							 		
	double  *layers_soilWaterConstant,
	double  *layers_soilWaterContent,
	double  *layers_extrSoilWaterCM,
	double  transpiration_proportion); // transpiration proportion = 1 for fallow

	
	
void soil_water(		
	int   	ndays, // number of days						
	double	*rainfall, // weather data
	double	*min_temp, 									
	double	*max_temp, 									
	double	*solar_radiation,						
	int     nlayers,  // soil parameters			
	double  upperlimitstage1evap,   	
	double	timeinstage2evap,			 		
	double	evapStage1, 						 		
	double	evapStage2, 	
	double	diffusConst, 						 		
	double	diffusSlope,	
	double  runoffCurveNumber,
	double  baresoilalbedo,					 		
	double	time, 
	double	tsw,
	double	tll,
	double	tsat,
	double  *layers_depth,
	double  *layers_airdry,	
	double  *layers_bd,
	double  *layers_wf,
	double  *layers_dryLowerLimit,
	double  *layers_drainedUpperLimit,
	double  *layers_saturation,							 		
	double  *layers_soilWaterConstant,
	double  *layers_soilWaterContent,
	double  *layers_extrSoilWaterCM,
	int     nstages, // number of crop stages
	int     *crop_stages_in_days, // index of the start date for crop ET stages,
	double  *crop_stages_coef,  // crop coeficient parameters
	double  *soilWaterContent); // return



#endif