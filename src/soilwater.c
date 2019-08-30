/********************************************************************

	soilwater.c

	Contains functions for soilwater modelling

	Author:         Fiona Evans
  Date:           September 2016
	Last updated    September 2016


*********************************************************************/
#include "include.h"

/*
 To be called by R function
*/

void soil_water_R(
    
    int 	  *ndays, // number of days						
    double	*rainfall, // weather data
    double	*min_temp, 									
    double	*max_temp, 									
    double	*radiation,						
    int     *nlayers,  // soil parameters	
    double  *upperlimitstage1evap,   	
    double	*timeinstage2evap,			 		
    double	*evapStage1, 						 		
    double	*evapStage2, 	
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
    double  *soilWaterContent)
{
  soil_water(*ndays, rainfall, min_temp, max_temp, radiation,
             *nlayers, *upperlimitstage1evap, *timeinstage2evap, *evapStage1, *evapStage2,
             *diffusConst, *diffusSlope, *runoffCurveNumber, *baresoilalbedo,
             *time, *tsw, *tll, *tsat,
             layers_depth, layers_airdry, layers_bd, layers_wf, 
             layers_dryLowerLimit, layers_drainedUpperLimit,
             layers_saturation, layers_soilWaterConstant, layers_soilWaterContent,
             layers_extrSoilWaterCM, *nstages, crop_stages_in_days, crop_stages_coef, soilWaterContent);
  
}





/* Calculate potential evaporation. */
double potential_evaporation(
    
  double min_temp,
  double max_temp,
  double solar_radiation,
  double baresoilalbedo)
{
  double weighted_temp, equilibriumET, potentialEvap;
  
  weighted_temp = 0.6 * max_temp + 0.4 * min_temp;
  equilibriumET = (solar_radiation / 0.04186) * 
									(0.000204 - 0.000183 * baresoilalbedo) * (weighted_temp + 29); 
  
  potentialEvap = equilibriumET * 1.1;    
  if (max_temp > 24) potentialEvap = equilibriumET * ((max_temp - 24) * 0.05 + 1.1) ;
  
  return potentialEvap;
}

/* 
Calculate actual evaporation for a given soil type using the potential
soil evaporation, upper limit of evaporation and amount of water infiltration.
Uses APSIM method hacked from APSIM Fortran code.
*/

double soilwat2_ritchie_evaporation(
    
  double  upperlimitstage1evap,   // soil parameters
  double	timeinstage2evap,			 
  double	*evapStage1, 						 
  double	*evapStage2, 						
  double	*time, 									 
  double	potentialEvap, 
  double	upperLimitEvap, 
  double	waterInfiltration)       // rainfall - runoff            
{
	
	double esoil1, esoil2, actualEvap;	
  
  //fprintf(stderr, "%f %f %f %f\n", waterInfiltration, *time, *evapStage1, *evapStage2);
	
	// if infiltration, reset evapStage1
  // reset evapStage2 if infil exceeds evapStage1
	if (waterInfiltration > 0) {
    *evapStage2 = MAX(0, *evapStage2 - MAX(0.0, waterInfiltration - *evapStage1)); 
    *evapStage1 = MAX(0, *evapStage1 - waterInfiltration); 
    
    // update time (in case evapStage2 changed)
    *time = pow(*evapStage2 / timeinstage2evap, 2);
  } 
  
  // are we in stage1 ?
  if (*evapStage1 < upperlimitstage1evap) {
    // we are in stage1, set esoil1 = potential, or limited by u.
    esoil1 = MIN(potentialEvap, upperlimitstage1evap - *evapStage1);
    
    if (potentialEvap > esoil1 && esoil1 < upperLimitEvap){    
      // potentialEvap not satisfied by 1st stage drying,
      // & there is evaporative sw excess to air_dry, allowing for esoil1.
      // need to calc. some stage 2 drying (esoil2).
      
      if (*evapStage2 > 0) {
        *time += 1.0;
        esoil2 = MIN(potentialEvap - esoil1, timeinstage2evap * sqrt(*time) - *evapStage2);
      }   
      else {				
				esoil2 = 0.6 * (potentialEvap - esoil1); // use ritchie's empirical transition constant (0.6)
			}
    }
    else{
      // no deficit (or esoil1.eq.eos_max,) no esoil2 on this day
      esoil2 = 0;
    }
    
    // check any esoil2 with lower limit of evaporative sw.
    esoil2 = MIN(esoil2, upperLimitEvap - esoil1);
    
    // update 1st and 2nd stage soil evaporation.  
    *evapStage1 += esoil1;
    *evapStage2 += esoil2;
    *time = pow(*evapStage2 / timeinstage2evap, 2);
    //fprintf(stderr, " In stage 1 %f %f ", *evapStage1, *evapStage2);
  }
  
  else{
    // no 1st stage drying. calc. 2nd stage    
    esoil1 = 0.0;    
    *time += 1;
    esoil2 = MIN(potentialEvap, timeinstage2evap * sqrt(*time) - *evapStage2);

    
    // check with lower limit of evaporative sw.    
    esoil2 = MIN(esoil2, upperLimitEvap);

    // update 2nd stage soil evaporation.    
    *evapStage2 += esoil2;
  }
  
  actualEvap = esoil1 + esoil2;
  
  // make sure we are within bounds
  if (actualEvap < 0) actualEvap = 0;
  if (actualEvap > potentialEvap) actualEvap = potentialEvap;
  if (actualEvap > upperLimitEvap) actualEvap = upperLimitEvap; 

  return actualEvap;
}

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
    double	waterInfiltration)       // rainfall - runoff            
{
  
  double esoil1, esoil2, se21, se22, actualEvap;	
  
  
  // if infiltration, reset evapStage1
  // reset evapStage2 if infil exceeds evapStage1
  if (waterInfiltration > 0) {
    *evapStage2 = MAX(0, *evapStage2 - MAX(0.0, waterInfiltration - *evapStage1)); 
    *evapStage1 = MAX(0, *evapStage1 - waterInfiltration); 
    
    // update time (in case evapStage2 changed)
    *time = pow(*evapStage2 / timeinstage2evap, 2);
  } 
  
  // are we in stage1 ?
  if (*evapStage1 < upperlimitstage1evap) {
    // we are in stage1, set esoil1 = potential, or limited by u.
    esoil1 = MIN(potentialEvap, upperlimitstage1evap - *evapStage1);
    esoil1 = MIN(0, MIN(esoil1, layers_soilWaterContent[0] - layers_airdry[0]));
    
    if (potentialEvap > esoil1){    
      // potentialEvap not satisfied by 1st stage drying,
      // & there is evaporative sw excess to air_dry, allowing for esoil1.
      // need to calc. some stage 2 drying (esoil2).
      
      if (*evapStage2 > 0) {
        *time += 1.0;
        esoil2 = MIN(potentialEvap - esoil1, timeinstage2evap * sqrt(*time) - *evapStage2);
      }   
      else {				
        esoil2 = 0.6 * (potentialEvap - esoil1); // use ritchie's empirical transition constant (0.6)
      }
      
      // Calculate stage two evaporation from layers 1 and 2
      se21 = MAX(0, MIN(esoil2, layers_soilWaterContent[0] - layers_airdry[0]));
      se22 = MAX(0, MIN(esoil2 - se21, layers_soilWaterContent[1] - layers_airdry[1]));
      
      esoil2 = se21 + se22;
    }
    else{
      // no deficit (or esoil1.eq.eos_max,) no esoil2 on this day
      esoil2 = 0;
    }
    
    
    // update 1st and 2nd stage soil evaporation.  
    *evapStage1 = upperlimitstage1evap;
    *evapStage2 += esoil2;
    *time = pow(*evapStage2 / timeinstage2evap, 2);
  }
  
  else{
    // no 1st stage drying. calc. 2nd stage and remove from laers 1 & 2  
    *evapStage1 = upperlimitstage1evap;
    *time += 1;
     
    esoil1 = 0.0;    
    esoil2 = MIN(potentialEvap, timeinstage2evap * sqrt(*time) - *evapStage2);
    
    se21 = MAX(0, MIN(esoil2, layers_soilWaterContent[0] - layers_airdry[0]));
    se22 = MAX(0, MIN(esoil2 - se21, layers_soilWaterContent[1] - layers_airdry[1]));
    
    esoil2 = se21 + se22;
    
    // update 2nd stage soil evaporation.    
    *evapStage2 += esoil2;
    
    //*evapStage2 = MAX(*evapStage2, 0);
    // fprintf(stderr, " In stage 2 %f %f \n", esoil2, *evapStage2);
  }
  
  actualEvap = esoil1 + esoil2;
  
  // make sure we are within bounds
  if (actualEvap < 0) actualEvap = 0;
  if (actualEvap > potentialEvap) actualEvap = potentialEvap;
 
  
  return actualEvap;
}


/*
Calculate runoff and remove from waterInfiltration.
*/
double runoff(

  int    nlayers,											// soil parameters
  double runoffCurveNumber,
  double tsw,
  double tsat,
  double *layers_wf,
  double *layers_soilWaterContent,
  double *layers_dryLowerLimit,
  double *layers_extrSoilWaterCM,
  double  waterInfiltration,         // weather data
  double  previousWaterInfiltration)
{
	int i;
	double runof, runof2, r2, winf1, twoDayRainfall, cn1, smx, sum, pb;
	
	// Was there rain?
  if (waterInfiltration > 0) {
     
    twoDayRainfall = waterInfiltration + previousWaterInfiltration;
    
    // Calculate runoff for large rainfall events
    if (twoDayRainfall >= 40) {    
      
      cn1 = -16.91 + (1.348 * runoffCurveNumber) -
        (0.01379 * pow(runoffCurveNumber,2)) +
        (0.0001172 * pow(runoffCurveNumber,3));
      smx = 254 * (100/cn1 - 1);
      
			sum = 0;
			for (i=0; i<nlayers; i++ ) {
				sum += layers_wf[i] * (layers_soilWaterContent[i] - layers_dryLowerLimit[i]) 
									/ layers_extrSoilWaterCM[i];
			}
			  
      r2 = smx * (1 - sum);
      if (r2 < 2.54) r2 = 2.54;
      
      pb = waterInfiltration - 0.2 * r2;
      runof = 0;
      if (pb > 0) runof = pow(pb,2) / (waterInfiltration + 0.8 * r2);
      
      winf1 = waterInfiltration - runof;
      
      if ((tsw * 10 + winf1) > (tsat * 10)) { 
        runof2 = tsw * 10 + winf1 - tsat * 10;
        runof = runof + runof2;   
        winf1 = winf1 - runof2;
      }
      waterInfiltration = winf1;
    }
    
  } 
  return waterInfiltration;
}

/*
Calculate downward flow and redistribute water through the profile.
*/
int flux(

  int    nlayers,
  double *layers_saturation,
  double *layers_soilWaterContent,
  double *layers_soilWaterConstant,
  double *layers_depth,
  double *layers_drainedUpperLimit,
  double  waterInfiltration)
{
	int idrsw = 0;
	double flux, hold, drainage;
	
	flux = waterInfiltration * 0.1; // Downward flux of water infiltration calculated in mm
  
	// Loop through layers
  for (int i=0; i<nlayers; i++) {
	
		hold = (layers_saturation[i] - layers_soilWaterContent[i]) * layers_depth[i];
		if (flux > hold) {
			drainage = layers_soilWaterConstant[i] * (layers_saturation[i] - layers_drainedUpperLimit[i]) *
									layers_depth[i];
			layers_soilWaterContent[i] = layers_saturation[i] - drainage / layers_depth[i];
			flux = flux - hold + drainage;
			idrsw = 1;
		}
		else {
			layers_soilWaterContent[i] = layers_soilWaterContent[i] + (flux / layers_depth[i]);
			if (layers_soilWaterContent[i] >= (layers_drainedUpperLimit[i] + 0.003)) {
				drainage = (layers_soilWaterContent[i] - layers_drainedUpperLimit[i]) * layers_soilWaterConstant[i] *
										layers_depth[i];
				layers_soilWaterContent[i] = layers_soilWaterContent[i] - (drainage / layers_depth[i]);
				flux = drainage;
				idrsw = 1;
			}
			else flux = 0;
		}	
	}
  drainage = flux * 10;
  
  return idrsw;
}


/*
Redistribute water between layers.  
*/
double flow(

  int idrsw,
  int nlayers,
  double diffusConst, 
  double diffusSlope,
  double *layers_soilWaterContent, 
  double *layers_dryLowerLimit, 
  double *layers_depth, 
  double *layers_saturation, 
  double *layers_drainedUpperLimit,
  double tll, 
  double *tsw) 
{
	double thet1, thet2,dbar, flow, swflow;
	double totExtrSoilWaterCM, soilWaterContentMM; 
	
	
	// Loop through layers
  for (int i=0; i<nlayers-1; i++) {

		thet1 = (layers_soilWaterContent[i] - layers_dryLowerLimit[i]); 
		thet2 = (layers_soilWaterContent[i+1] - layers_dryLowerLimit[i+1]);  
		if (thet1 < 0) thet1 = 0 ;
		if (thet2 < 0) thet2 = 0;
		
		dbar = diffusConst/100 * exp(diffusSlope * (thet1 + thet2)/2); // cm/day
		if (dbar > 100) dbar = 100;
		
		// flow = amount of upward flow between two layers = dbar * volumetric soil water gradient
		flow = dbar * (thet2 - thet1)/((layers_depth[i] + layers_depth[i+1])/2);
		
		swflow = flow / layers_depth[i];
		
		if (flow > 0) {
			if (idrsw == 1)  {
			  if ((layers_soilWaterContent[i] + swflow) > layers_saturation[i]) {
				  flow = (layers_saturation[i] - layers_soilWaterContent[i]) * layers_depth[i];
			  }
			}
			if (idrsw == 0) {
			  if ((layers_soilWaterContent[i] + swflow) > layers_drainedUpperLimit[i]) {
				  flow = (layers_drainedUpperLimit[i] - layers_soilWaterContent[i]) * layers_depth[i];
			  }
			}
		}
		if (flow <= 0) {
			if (idrsw == 1) {
			  if  ((layers_soilWaterContent[i+1] - swflow) > layers_saturation[i+1]) {
				  flow = (layers_soilWaterContent[i+1] - layers_saturation[i+1]) * layers_depth[i+1];
			  }
			}
				
			if (idrsw == 0) {
			  if ((layers_soilWaterContent[i+1] - swflow) > layers_drainedUpperLimit[i+1]){
				  flow = (layers_soilWaterContent[i+1] - layers_drainedUpperLimit[i+1]) * layers_depth[i+1];
			  }			
			} 
		
		} 
		
		layers_soilWaterContent[i] += (flow / layers_depth[i]);
		layers_soilWaterContent[i+1] -= (flow / layers_depth[i+1]) ;
		
		layers_soilWaterContent[i] = MAX(layers_soilWaterContent[i], 0); // cm
		layers_soilWaterContent[i+1] = MAX(layers_soilWaterContent[i+1], 0); // cm  
		}
  

	*tsw = 0;
	for (int i=0; i<nlayers; i++) {
		*tsw += layers_soilWaterContent[i] * layers_depth[i];
	  }
	if (*tsw < tll) *tsw = tll;
	 
	totExtrSoilWaterCM = *tsw - tll;
  soilWaterContentMM = totExtrSoilWaterCM * 10;  
  
  return soilWaterContentMM;
}



/* 
Iterate soil water model for one day. 

The water balance calculation method is derived from Ritchie, J. T. 
(1972). Model for predicting evaporation from a crop with incomplete 
cover. Water Resour. Res., 8, 1204?1213

FAO calculation of crop evapotranspiration from
'http://www.fao.org/docrep/x0490e/x0490e0b.htm#chapter 6'
*/
double soil_water_day_evapotranspiration( 
																							
  double	waterInfiltration,          // weather data
  double	previousWaterInfiltration,  
  double	min_temp, 									
  double	max_temp, 									
  double	solar_radiation,						
  int     nlayers,                    // number of soil layers																							
  double  upperlimitstage1evap,   		// soil parameters
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
  double  transpiration_proportion)
{
	int idrsw = 0;
	double upperLimitEvap, actualEvap, potentialEvap, soilWaterContentMM, plantET, tmp; 
	
	// Calculate and remove runoff
	waterInfiltration = runoff(nlayers, runoffCurveNumber, *tsw, tsat, layers_wf,
                            layers_soilWaterContent, layers_dryLowerLimit, layers_extrSoilWaterCM,
                            waterInfiltration, previousWaterInfiltration);
	
	// Distribute downward flux through layers
	idrsw = flux(nlayers, layers_saturation, layers_soilWaterContent, layers_soilWaterConstant,
              layers_depth, layers_drainedUpperLimit, waterInfiltration);
  
  
  // Potential evaporation
  potentialEvap = potential_evaporation(min_temp, max_temp, solar_radiation, 
						baresoilalbedo);
  
						
	// Plant evapotranspiration
	plantET = transpiration_proportion * potentialEvap;
	
	// Actual evaporation 

	upperLimitEvap = MAX(0, (layers_soilWaterContent[0] - layers_airdry[0]) * layers_depth[0]); // mm/day
	upperLimitEvap += MAX(0, (layers_soilWaterContent[1] - layers_airdry[1]) * layers_depth[1]);
	if (upperLimitEvap > potentialEvap) upperLimitEvap = potentialEvap;
	
	actualEvap = soilwat2_ritchie_evaporation(upperlimitstage1evap, timeinstage2evap, evapStage1, 
                                           evapStage2, time, potentialEvap, 
                                           upperLimitEvap, waterInfiltration); 
	
	
	
	
	

  /*
  actualEvap = soilwat_howleaky_evaporation(upperlimitstage1evap, timeinstage2evap, evapStage1, 
                                            evapStage2, time, potentialEvap, 
                                            nlayers, layers_soilWaterContent, layers_airdry,
                                             waterInfiltration); 
  */
   
  
  // Problem with evapStage2 when there is plantET ?
	
	actualEvap = plantET + actualEvap;
	
	//fprintf(stderr, "%f %f %f\n", layers_soilWaterContent[0], layers_soilWaterContent[1], actualEvap * 0.1 / layers_depth[0]);
	
	// Extract evaporation from layer 1
	if (actualEvap * 0.1 / layers_depth[0] < layers_soilWaterContent[0]) {
	  layers_soilWaterContent[0] = layers_soilWaterContent[0] - actualEvap  * 0.1 / layers_depth[0];
	}
	else {
	  layers_soilWaterContent[0] = 0;
	  tmp = (actualEvap * 0.1 / layers_depth[0]) - layers_soilWaterContent[0];
	  tmp = tmp * layers_depth[0] / 0.1;
	  layers_soilWaterContent[1] = layers_soilWaterContent[1] - tmp  * 0.1 / layers_depth[1];
	}
	  
	
  
  // Redistribute water between layers  (flow = amount of upward flow between two layers)
  soilWaterContentMM = flow(idrsw, nlayers, diffusConst, diffusSlope,
					   layers_soilWaterContent, layers_dryLowerLimit, layers_depth, 
						 layers_saturation, layers_drainedUpperLimit,  tll, tsw); 
  
  
  return soilWaterContentMM;
}

/* Run fallow soil water model for a sequence of days */
void soil_water( 
				
  int 	  ndays, // number of days																							
  double	*rainfall, // weather data
  double	*min_temp, 									
  double	*max_temp, 									
  double	*radiation,						
  int     nlayers,  // number of soil layers																							
  double  upperlimitstage1evap,   		// soil parameters
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
  double  *soilWaterContent) // return
{
	double	waterInfiltration;          // weather data
	double	previousWaterInfiltration;
	double  transpiration_proportion = 0;
	
	time = 0;
	
	waterInfiltration = rainfall[0];
	previousWaterInfiltration = 0;
	
	soilWaterContent[0] = soil_water_day_evapotranspiration(waterInfiltration, previousWaterInfiltration, 
			min_temp[0], max_temp[0], radiation[0],
			nlayers, upperlimitstage1evap, timeinstage2evap, &evapStage1, &evapStage2,
			diffusConst, diffusSlope, runoffCurveNumber, baresoilalbedo,
			&time, &tsw, tll, tsat,
			layers_depth, layers_airdry, layers_bd, layers_wf, layers_dryLowerLimit, layers_drainedUpperLimit,
			layers_saturation, layers_soilWaterConstant, layers_soilWaterContent,
			layers_extrSoilWaterCM, transpiration_proportion=1);
	
	for (int i=1; i<ndays; i++){
	  
	  transpiration_proportion = 0;
	  for (int j=0; j<nstages; j++) {
	    if (i >= crop_stages_in_days[j]) {
	      transpiration_proportion=crop_stages_coef[j];
	    }
	  }
	  
	 // fprintf(stderr, "%f\n", transpiration_proportion);
	  
		waterInfiltration = rainfall[i];
		previousWaterInfiltration = rainfall[i-1];
		
		soilWaterContent[i] = soil_water_day_evapotranspiration(waterInfiltration, previousWaterInfiltration, 
				min_temp[i], max_temp[i], radiation[i],
        nlayers, upperlimitstage1evap, timeinstage2evap, &evapStage1, &evapStage2,
        diffusConst, diffusSlope, runoffCurveNumber, baresoilalbedo,
        &time, &tsw, tll, tsat,
        layers_depth, layers_airdry, layers_bd, layers_wf, layers_dryLowerLimit, layers_drainedUpperLimit,
        layers_saturation, layers_soilWaterConstant, layers_soilWaterContent,
        layers_extrSoilWaterCM, transpiration_proportion);
		
	}
	
	
	
}

