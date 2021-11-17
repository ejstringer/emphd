
# fst within trips -------------------------------------------------------------

#' Within trip fst calculations
#' 
#' fst bootstraps df calculated for a single genlight object.
#' fst separated by trip and based on your defined populations.
#' I have also included minimum sample size. You no longer need to
#' 
#' 
#' @param glx -- a genlight object
#' @param as.pop -- defined population to base fst calculations on
#' @param min.n -- minimum sample size
#' @param equal -- equal sample size
#' @return bootstraps df from gl.fst.pop
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples fst <- em.fst.trip(gl, gridId, 4)

em.fst.trip <- function(glx, as.pop = "gridId", min.n = 4, equal = F){
  
      
      # subset for fst calculations
      tripPop <- paste0(glx@other$ind.metrics$trip, 
                        glx@other$ind.metrics[,as.pop])
      
      glTidyn <- em.gl.min(glx, define.pop = tripPop, s.size = min.n, 
                equal = equal)
      
      glTidy <- em.gl.remove.singlepop(glTidyn,define.pop = as.pop)
      
      
      # seperate by trip
      pop(glTidy) <- glTidy@other$ind.metrics$trip
      glseptrip <- seppop(glTidy)
      
      glseptrip  <- lapply(glseptrip, 
                           function(x) em.gl.assignPop(x, define.pop = as.pop))
      
      
      # calculate fst
      fst <- pbapply::pblapply(glseptrip,
                               function(x) gl.fst.pop(x, ncluster=4))
      
      boot.fst <- lapply(fst, function(x) x$Bootstraps)
      
      for (i in seq(length(boot.fst))) {
        boot.fst[[i]]$trip <- names(boot.fst)[[i]]
      }
      
      
      # sample size
      sampleSize <- data.frame(table(glTidy$other$ind.metrics[, as.pop], 
                          glTidy$other$ind.metrics$trip))
      sampleSize2 <- sampleSize 
      names(sampleSize) <- c("Population1", "trip", "n1")
      names(sampleSize2) <- c("Population2", "trip", "n2")

      
      fst.reduce <-  do.call("rbind", boot.fst)
      fst.n <- merge(fst.reduce, sampleSize, by = c("trip", "Population1"))
      fst.nn <- merge(fst.n, sampleSize2, by = c("trip", "Population2"))
      
      fst <- fst.nn
      names(fst) <- c(names(fst.nn)[1:3], paste0(1:100, "boot"), 
                      names(fst.nn)[104:length(fst.nn)])
      fstboots <- fst[,paste0(1:100, "boot")]
      fstvar <- apply(fstboots, MARGIN = 1, var)
      fst.nn$weights <- 1/fstvar

  return(fst.nn)
}

# fis within trips -------------------------------------------------------------

#' Within trip fis calculations
#' 
#' fis bootstraps df calculated for a single genlight object.
#' fis separated by trip and based on your defined populations.
#' 
#' @param glx -- a genlight object
#' @param as.pop -- defined population to base fis calculations on
#' @param min.n -- minimum sample size
#' @param equal -- equal sample size
#' @return df with he, ho, and fis
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples  heterozygosity <- em.fis.trip(gl)

em.fis.trip <- function(glx, as.pop = "gridId", min.n = 4, equal = F){
  
  # subset by sample size
  tripPop <- paste0(glx@other$ind.metrics$trip, 
                    glx@other$ind.metrics[,as.pop])
  
  glTidyn <- em.gl.min(glx, define.pop = tripPop, s.size = min.n, 
                       equal = equal)
  if(is.null(glTidyn)) stop("min.n is larger than max number of samples")
  # sample size
  sampleSize <- data.frame(table(glTidyn$other$ind.metrics[, as.pop], 
                                 glTidyn$other$ind.metrics$trip))
  names(sampleSize) <- c(as.pop, "trip", "n")
   

  
  
  # seperate by trip and as.pop
  pop(glTidyn) <- paste(glTidyn@other$ind.metrics$trip, 
                       glTidyn@other$ind.metrics[,as.pop])
  glsep <- seppop(glTidyn)
  glsepHo <- sapply(lapply(glsep, gl.Ho), function(x) mean(x, na.rm = T))
  glsepHe <- sapply(lapply(glsep, gl.He), function(x) mean(x, na.rm = T))
  glsepHoVar <- sapply(lapply(glsep, gl.Ho), function(x) var(x, na.rm = T))
  glsepHeVar <- sapply(lapply(glsep, gl.He), function(x) var(x, na.rm = T))
  
  dfnames <- data.frame(do.call("rbind", strsplit(names(glsep), "\\s+")),
                        row.names = NULL)
  names(dfnames) <- c("trip", as.pop)
  df.n <- merge(dfnames, sampleSize, by = c("trip", as.pop))
  
  heterozygosity <- data.frame(df.n, row.names = NULL,
                               ho = glsepHo, hoVar = glsepHoVar,
                               he = glsepHe, heVar = glsepHeVar,
                               fis = (glsepHe- glsepHo)/glsepHe)
  return(heterozygosity)
}

