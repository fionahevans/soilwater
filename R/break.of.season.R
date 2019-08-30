

#' Run the break of season
#' 
#' Follows a simple 2-pert rule: break occurs if there is 25mm of rainfall over three 
#' days after 25 April, or if there is 5mm over 3 days after 5 June 
#'
#' @param r Data frame containing daily climate data.
#' @param year Year we are considering.
#'
#' @return Date of break of season.
#' @export
break.of.season <- function(r, year, mm1=25, mm2=5, 
                            date1=as.Date(paste("2504", year, sep=""), "%d%m%Y"),
                            date2=as.Date(paste("0506", year, sep=""), "%d%m%Y")) { 
  
  
  if (r[nrow(r), "DATE"] < as.Date(paste("2504", year, sep=""), "%d%m%Y")) {
    ret <- r[nrow(r), "DATE"]
  }
  
  else {
    
    r <- subset(r, DATE <= as.Date(paste("3110", year, sep=""), "%d%m%Y"))
    
    # 25mm over three days after 25 April, then 5mm over 3 days after 5 Jun
    r[, "rsum"] <- sumFilter(r[, "RAIN"], 3)
    
    sub <- subset(r, DATE >= as.Date(paste("2504", year, sep=""), "%d%m%Y") &
                    DATE < as.Date(paste("0506", year, sep=""), "%d%m%Y"))
    
    j <- which(sub[, "rsum"] > mm1)[1]
    
    if (is.na(j)){
      sub <- subset(r, DATE >= date2 )
      j <- which(sub[, "rsum"] > mm2)[1]
    }
    
    if (is.na(j)) j <- nrow(sub)
    ret <- sub[j, "DATE"]
  }
  ret
}

sumFilter <- function(x, n) { 
  # Sum n elements of vector x and store sum of last three days
  
  ret <- rep(NA, length(x))  
  
  for (i in n:length(x)){
    ret[i] <- sum(x[(i-n+1):i])  
  }
  
  round(ret, 1)
}