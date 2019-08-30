

#' Run fallow soil water model
#' 
#' The water balance calculation method is derived from Ritchie, J. T. 
#' (1972). Model for predicting evaporation from a crop with incomplete 
#' cover. Water Resour. Res., 8, 1204?1213
#'
#' @param climate.data Data frame containing daily climate data.
#' @param soil.init Data frame containing soil data.
#'
#' @return Vector of daily plant available soil water.
#' @useDynLib soilwater
#' @export
soil.water <- function(climate.data, soil.init, season.break=NULL, wheat.et.table=NULL){

	ndays <- nrow(climate.data)
	rainfall <- climate.data[, "RAIN"]
	min_temp <- climate.data[, "MIN_TEMP"]
	max_temp <- climate.data[, "MAX_TEMP"]
	radiation <- climate.data[, "RADIATION"]
	
	nlayers <- 2
	
	# Fallow model
	nstages <- 0
	stages.days <- rep(0, 3)
	stages.coef <- rep(1, 3)

	
	# Crop model
	if (!is.null(wheat.et.table) && !is.null(season.break)) {
	  nstages <- nrow(wheat.et.table)
	  stages.days <- rep(season.break, nstages)
	  for (j in 1:(nstages-1))  {
	    stages.days[j+1] <- sum(wheat.et.table[1:j, "length"]) + season.break
	  }
	  stages.coef <- wheat.et.table[, "coef"]
	}
	
	#cat("stages.days = ", stages.days, "\n");
	#cat("stages.coef = ", stages.coef, "\n");

	
  soilWaterContentMM <- .C("soil_water_R", 	
			as.integer(ndays),																							
			as.double(rainfall),
			as.double(min_temp), 									
			as.double(max_temp), 									
			as.double(radiation),						
			as.integer(nlayers),																					
			as.double(soil.init["UPPERLIMITSTAGE1EVAP"]),
			as.double(soil.init["TIMEINSTAGE2EVAP"]),			 		
			as.double(soil.init["EVAPSTAGE1"]), 						 		
			as.double(soil.init["EVAPSTAGE2"]),	
			as.double(soil.init["DIFFUSCONST"]), 						 		
			as.double(soil.init["DIFFUSSLOPE"]), 	
			as.double(soil.init["RUNOFFCURVENUMBER"]),
			as.double(soil.init["BARESOILALBEDO"]),				 		
			as.double(soil.init["TIME"]), 
			as.double(soil.init["TSW"]),
			as.double(soil.init["TLL"]),
			as.double(soil.init["TSAT"]),
			as.double(c(soil.init["LAYER1_DEPTH"], soil.init["LAYER2_DEPTH"])),
			as.double(c(soil.init["LAYER1_AIRDRY"], soil.init["LAYER2_AIRDRY"])),	
			as.double(c(soil.init["LAYER1_BD"], soil.init["LAYER2_BD"])),
			as.double(c(soil.init["LAYER1_WF"], soil.init["LAYER2_WF"])),
			as.double(c(soil.init["LAYER1_DRYLOWERLIMIT"], soil.init["LAYER2_DRYLOWERLIMIT"])),
			as.double(c(soil.init["LAYER1_DRAINEDUPPERLIMIT"], soil.init["LAYER2_DRAINEDUPPERLIMIT"])),
			as.double(c(soil.init["LAYER1_SATURATION"], soil.init["LAYER2_SATURATION"])),								 		
			as.double(c(soil.init["SOILWATERCONSTANT"], soil.init["LAYER2_SOILWATERCONSTANT"])),
			as.double(c(soil.init["LAYER1_SOILWATERCONTENT"], soil.init["LAYER2_SOILWATERCONTENT"])),
			as.double(c(soil.init["LAYER1_EXTRSOILWATERCM"], soil.init["LAYER2_EXTRSOILWATERCM"])),
			as.integer(nstages),
			as.integer(stages.days),
			as.double(stages.coef),
			soilWaterContent = double(ndays),
			PACKAGE="soilwater")$soilWaterContent
			
  soilWaterContentMM		
}
