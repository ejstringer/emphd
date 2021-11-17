
# grid or site distance --------------------------------------------------------

#' distance between grids or sites
#' 
#' simple function to create a dataframe of distance in metres 
#' between sites
#' 
#' @param dfcoor -- a df with lat and lon, and location id
#' @param between -- the location id to calculate between distances
#' @param proj -- projection of your lat lons
#' @param transform -- projection to transform into (UTM)
#' @return df of distances between 'populations'
#' @export
#' @author Emily Stringer
#' @examples em.distance.m(grids, between = "gridId")

em.distance.m <- function(dfcoor, between = "site", proj = "+proj=longlat",
                        transform = "+init=epsg:32754"){
  cord.dec = sp::SpatialPoints(cbind(dfcoor$lon, dfcoor$lat), 
                               proj4string= sp::CRS(proj))
  cord.UTM <- sp::spTransform(cord.dec, sp::CRS(transform))  
  
  mat <- as.matrix(dist(cord.UTM@coords))
  matNames <- dfcoor[, between]
  rownames(mat) <- matNames
  colnames(mat) <- matNames
  
  distmat <- as.dist(mat)
  
  dfdist <- otuSummary::matrixConvert(distmat, 
                                      colname = c("Population1", 
                                                  "Population2", 
                                                  "metres")) 

  return(dfdist)
  
}
