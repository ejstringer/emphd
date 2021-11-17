# assign pop -------------------------------------------------------------------

#' reassigne population 
#' 
#' I created this function to help me with reassigning pop()
#' for a list of gl objects, using lapply.
#' 
#' 
#' @param glx -- genlight object
#' @param define.pop -- "variable" from ind.metrics to assign as new pop
#' @return genlight object with new pop() values
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples 
#' phyear <- seppop(gl) %>% 
#' lapply(., function(x) em.gl.assignPop(x, define.pop = "pop"))


em.gl.assignPop <- function(glx, define.pop = "gridId"){
  
  if (length(define.pop) == 1){
    if(!define.pop %in% names(glx@other$ind.metrics)){
      stop(paste(define.pop, 'not in ind.metrics; names(glx@other$ind.metrics)')) 
      }
    pop(glx) <- glx@other$ind.metrics[,define.pop]
  }else{  
    pop(glx) <- define.pop
  }
  
  return(glx)
}

# min sample per pop -----------------------------------------------------------

#' minimum population size
#' 
#' removes populations that have less than defined sample size
#' 
#' 
#' @param glx -- genlight object
#' @param define.pop -- "variable" from ind.metrics to assign as pop
#' @param s.size -- minimum sample size
#' @param equal -- if each pop should have the same number of samples
#' @return genlight object
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples 
#' em.gl.min(pherm, "gridId", 4, FALSE)

em.gl.min <- function(glx, define.pop = "gridId", s.size = 4, equal = FALSE){
  
  
  my.fun <- function(x, equal){
    if (length(x@ind.names) >= s.size) {
      aa <- x@ind.names
      if(equal) aa <- sample(x$ind.names, s.size, replace = FALSE)

    }else{
      aa <- NULL
    }
    return(aa)
  }

  glx <- em.gl.assignPop(glx, define.pop)
  glsep <- seppop(glx)
  
  indexNames <-  unlist(lapply(glsep, function(x) my.fun(x, equal)))
  
  glSample <- NULL
  if(length(indexNames) > 0) glSample <- gl.keep.ind(glx, ind.list = indexNames)
    
  return(glSample)
}


# remove trips for fst ---------------------------------------------------------

#' remove single population trips
#' 
#' removes trips that only have one pop(), where fst's can't be computed
#' 
#' 
#' @param glx -- genlight object
#' @param define.pop -- "variable" from ind.metrics to assign as pop
#' @return genlight object with removed trips
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples 
#' phMinSample <- lapply(1:5, function(x) em.gl.min.pop.n(pherm, 
#'                                 define.pop = iNewPop, s.size = x))

em.gl.remove.singlepop <- function(glx, define.pop = "gridId"){
  if(!define.pop %in% names(glx@other$ind.metrics)){
    stop(paste(define.pop, 'not in ind.metrics; names(glx@other$ind.metrics)')) 
  }
  if(!"trip" %in% names(glx@other$ind.metrics)) stop('trip not in ind.metrics') 
  
  
  npopsTrip <- rowSums(table(glx@other$ind.metrics$trip, 
                             glx@other$ind.metrics[,define.pop]) != 0)
  goodTrip <- names(npopsTrip)[npopsTrip > 1]
  
  gl2pops <- gl.keep.pop(glx, pop.list = goodTrip, as.pop = "trip")  
  return(gl2pops)
  
}

