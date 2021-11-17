# filter -----------------------------------------------------------------------

#' My Filtering protocol
#' 
#' Currently I am filtering on call rates, monomorphs, 
#' minor allele frequency, and secondaries.
#' 
#' @param glx -- a genlight object
#' @param t.Loc -- callrate threshold for loci
#' @param t.Ind -- callrate threshold for individuals
#' @param t.maf -- threshold for minor allele frequency
#' @return my filtered genlight object
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples gl <- em.filtering.gl(gl)

em.filtering.gl <- function(glx, t.Loc = 0.95, t.Ind = 0.95, t.maf = 0.005){
  
  glO <- glx
  
    glx <-  gl.filter.callrate(glx ,method = "loc", threshold = t.Loc)
    
    glx <- gl.filter.callrate(glx, method="ind", threshold = t.Ind)
    
    glx <- gl.filter.monomorphs(glx)
    
    glx <- gl.filter.maf(glx,threshold = t.maf)
    
    glx <- gl.filter.secondaries(glx) #?? should I filter for secondaries
 
cat("\n"); cat("\n")  
cat("pre filtering loc:", nLoc(glO), "and ind:",  nInd(glO), "\n")
cat("post filtering loc:", nLoc(glx), "and ind:",  nInd(glx), "\n")

return(glx)

}
# history ----------------------------------------------------------------------

#' My Filtering history
#' 
#' Shows how I filtered my data in a nice df (usable in kable and DT)
#' 
#' 
#' @param glx -- a filtered genlight object
#' @return filtering history as a df
#' @export
#' @import dartR
#' @author Emily Stringer
#' @examples em.filtering.history(gl)

em.filtering.history <- function(glx){
  history <- data.frame(nr = 1:length(as.character(glx@other$history)), 
                       history = as.character(glx@other$history)) 
  
  history$history <- as.character(history$history)
  history <- rbind(history, c("->", as.character(substitute(glx))))
    
  return(history)
}



