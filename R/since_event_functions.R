# peaks/mixing trips -----------------------------------------------------------

#' Time sincePeak
#' 
#' This function takes dates and sets the time since peak to zero at 
#' those dates and then every month after is +1. It can now be used on
#' any dataframe with a 'trip' column and date format yyyy-mm-dd 
#' 
#' @param eventdate -- date of event
#' @param df -- df with variable trip
#' @return df with updated sincePeak and peakNo
#' @export
#' @author Emily Stringer 
#' @examples 
#' em.months.since(c("2004-02-01", "2010-01-01"), df)

em.months.since <- function(eventdate, df){
  if(!"trip" %in% names(df)) stop('trip not in ind.metrics') 
  
  eventdate <- sort(eventdate)
  mintrip <- min(lubridate::ymd(df$trip))
  maxtrip <- max(lubridate::ymd(df$trip))
  
  dfdates <- data.frame(trip=seq(mintrip, maxtrip, by="months"), 
                        monthsSince = NA, sinceEvent = NA)
  
  
  for (i in 1:length(eventdate)) {
    p2 <- grep(eventdate[i],dfdates$trip) # date of event (zero)
    
    if(length(p2) > 0){
      sp <- (1:(nrow(dfdates)-p2+1))-1 # since event
      dfdates[p2:nrow(dfdates),]$monthsSince <- sp
      dfdates[p2:nrow(dfdates),]$sinceEvent <- rep(paste0("event", i), 
                                                   length(sp))
      
    }
  }
  
  df$trip <- as.character(df$trip)
  dfdates$trip <- as.character(dfdates$trip)
  
  
  dfsince <- merge(df, dfdates, by = "trip")
  
  return(dfsince)
}

