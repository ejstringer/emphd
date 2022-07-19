
# boom bust sim ----------------------------------------------------------------
#' conceptual simulation gl.offspring 
#' 
#' This function creates a simulation using gl.sim.offspring
#' over 20 generations, where five populations are isolated 
#' except during a few generations when dispersal happens. It 
#' also simulates population increases during across two 
#' generations, the first population increase with dispersal, 
#' the second the next generation after dispersal. I plan on 
#' developing a simulation with constant dispersal as well.
#' 
#' Make sure sex is defined by m and f (only these will be included) 
#' 
#' @param glBase -- a genlight object to base simulation on (starting snps)
#' @param pops -- generations when mixing occurs
#' @param popSize -- number of individuals per population during bust
#' @param dispersal -- amount of dispersal
#' @return returns genlight for 20 generations
#' @export
#' @import tidyverse
#' @import dartR
#' @author Emily Stringer
#' @examples 
#' em.sim.boombust(pherm, pops = 5, popSize = 50, dispersal = 0.25, 
#'                                                  dispConstant = TRUE)



em.sim.boombust <- function(glBase, pops = 5, popSize = 50, 
                            dispersal = 0.25, dispConstant = FALSE) {
  
  #tidy sex column 
  if(is.null(glBase@other$ind.metrics$sex)) stop("base gl needs sex column")
  sex <- glBase@other$ind.metrics$sex
  sexId <- unique(sex[c(grep("m", sex, ignore.case = T),
  grep("f", sex, ignore.case = T))])
  
  
  glsex <- gl.keep.pop(glBase, pop.list = sexId, as.pop = "sex")
  sex <- glsex@other$ind.metrics$sex
  sexTidy <- factor(ifelse(grepl("f", sex, ignore.case = T), "f", "m"))
  glsex@other$ind.metrics$sex <- sexTidy
  
  glfilter <- gl.filter.callrate(glsex, threshold = 1, plot = F) # filter NAs
  
  #number of populations
  n.pop= pops
  popSize = popSize
  #number of generations to run
  n.gen= 20
  
  m.rate = 1e-7 #mutation rate
  mixev = c(5,15,21) #mix generatoion
  if (dispConstant) mixev = 1:n.gen
  pop.inc = c(5,15,21) 
  m.ind = dispersal #mixing ratio - just going with 1 ind
  calcev = 1 #calc fst every calcev generation DONT NEED
  
  #start with a list of population based on pstart/glBase
  p <- list()
  pstart <- glfilter
  
  
  for (j in 1:n.pop){
    dads <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "m")
    mums <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "f")
    
    p[[j]] <- gl.sim.offspring(fathers = dads, mothers = mums,
                               noffpermother = 4, sexratio = 0.5)
    n.ind <- nInd(p[[j]])
    p[[j]]@pop <- factor(rep(paste0("P", j), n.ind))
    
    indNames(p[[j]]) <- paste0("p",j,"_",1:n.ind, "_",sprintf("%02d", 0))
  } #end for loop
  
  
  #generation zero
  em.gl.indmetrics <- function(glsim, gen = 0, disp.rate = 0){
    
    genn <- paste0("gen", sprintf("%02d", gen))
    glsim@other$ind.metrics <- data.frame(id = glsim@ind.names,
                                          pop = glsim@pop,
                                          trip = genn,
                                          dispersal = disp.rate, 
                                          sex = glsim@other$sex)
    return(glsim@other$ind.metrics)
    
  }
 
  glList <- list()
  pIndMet <- do.call("rbind", lapply(p, em.gl.indmetrics))
  pp <- do.call("rbind", p)
  pp@other$ind.metrics <- pIndMet
  
  cc <- 1
  glList[[cc]] <- pp 
  
  
  #number of mixing combinations between pops
  pcomb <- t(combn(n.pop,2))
  
  
  for(i in 1:n.gen){
    
    
    if ((i %in% mixev)){ #mix every mixev generations
      
      for (iii in 1:nrow(pcomb)){
        pa <- pcomb[iii,1]
        pb <- pcomb[iii,2]
        
        
        n.inda <- nInd(p[[pa]]) # different number of individuals
        n.indb <- nInd(p[[pb]]) # per population
        
        nmi <- min(n.indb, n.inda)*m.ind #number of individuals that swap
        
        pam <- sample(1:n.inda, nmi, replace = F)
        pbm <- sample(1:n.indb, nmi, replace = F)
        pback <- rbind(p[[pa]][-pam,], p[[pb]][pbm,])
        
        
        sexA <- factor(c(as.character(p[[pa]][-pam,]@other$sex),
                         as.character(p[[pb]][pbm,]@other$sex))) 
        pback@other$sex <- sexA                          
        popNames(pback) <- c(paste0("P",pa),paste0("P",pa))
        
        
        
        pback2 <- rbind(p[[pb]][-pbm,], p[[pa]][pam,]) 
        
        sexB <- factor(c(as.character(p[[pb]][-pbm,]@other$sex), 
                         as.character(p[[pa]][pam,]@other$sex)))
        pback2@other$sex <- sexB
        
        
        popNames(pback2) <- c(paste0("P",pb),paste0("P",pb))
        p[[pa]] <- pback
        p[[pb]] <- pback2
      }
    }
    
    
    for (ii in 1:n.pop){ # they breed creating the next generation
      
      gen <- p[[ii]]
      
      genDad <- gen[gen@other$sex == "male",]
      
      genMum <- gen[gen@other$sex == "female",]
      
      gen <- gl.sim.offspring(fathers = genDad, mothers = genMum,
                              noffpermother = 4, sexratio = 0.5)
      
      genF <- gen[gen@other$sex == "female",]
      genF25 <- genF[sample(1:nInd(genF), min(nInd(genF), popSize/2)),]
      if ((i %in% pop.inc | (i-1) %in% pop.inc)) genF25 <- genF
      
      genM <- gen[gen@other$sex == "male",]
      genM25 <- genM[sample(1:nInd(genM), min(nInd(genM), popSize/2)),]
      if ((i %in% pop.inc | (i-1) %in% pop.inc)) genM25 <- genM
      
      
      gen50 <- rbind(genM25, genF25)
      gen50@other$sex <- factor(c(as.character(genM25$other$sex),
                                  as.character(genF25$other$sex)))
      
      
      n.ind <- nInd(gen50)
      gen50@pop <- factor(rep(paste0("P", ii), n.ind))
      indNames(gen50) <- paste0("p",ii,"_",
                                sprintf("%02d", 1:n.ind),
                                "_",sprintf("%02d", i))
      
      p[[ii]] <- gen50
      
      
    }
    
    # save generation
      
      
      dr <- ifelse(i %in% mixev, dispersal, 0)
      p2 <- lapply(p, function(x) em.gl.indmetrics(x, gen = i, disp.rate = dr))
      pIndMet <- do.call("rbind", p2)
      pp <- do.call("rbind", p)
      pp@other$ind.metrics <- pIndMet
      
      cc <- cc+1
      glList[[cc]] <- pp 
      
      
      
      cat(paste("Generation:", i,"/" ,n.gen, "\n"))
      flush.console()
  
    
  }
   
  
   indmet <- do.call("rbind", lapply(glList, function(x) x@other$ind.metrics))
   res <- do.call("rbind", glList)
   res@other$ind.metrics <- indmet
   res@other$loc.metrics.flags <- glBase@other$loc.metrics.flags
   
   
  return(res)
  
}


# density dependent dispersal --------------------------------------------------


#' conceptual simulation gl.offspring in regard
#' 
#' quick and dirty sim function based on em.sim.boombust.
#' looking at density-dependent dispersal.
#' 
#' d-d = density dependent.
#' See em.sim.boombust for more details...
#' 
#' 
#' @param glBase -- a genlight object to base simulation on (starting snps)
#' @param popSize -- subpopulaiton size
#' @param subpops -- number of subpopulations
#' @param dispersal -- amount of dispersal
#' @param dispersalType -- positive, negative, or constant dispersal.
#' @param nOff -- number of offspring after generation 0 (gen0 is set to 4)
#' @return returns genlight for 20 generations
#' @export
#' @import tidyverse
#' @import dartR
#' @author Emily Stringer
#' @examples sim <- em.sim.boombust.ddd(gl)


em.sim.boombust.ddd <- function(glBase, popSize = 20, subpops = 10, 
                                dispersal = 0.10, 
                                dispersalType = "positive", nOff = 4) {
  
  dispOpt <-  c("positive", "negative", "constant")
  if(!dispersalType %in%  dispOpt) stop("dispersalType must equal either:", 
                                       paste(dispOpt, collapse = " "))
  #tidy sex column 
  if(is.null(glBase@other$ind.metrics$sex)) stop("base gl needs sex column")
  sex <- glBase@other$ind.metrics$sex
  sexId <- unique(sex[c(grep("m", sex, ignore.case = T),
                        grep("f", sex, ignore.case = T))])
  
  
  glsex <- gl.keep.pop(glBase, pop.list = sexId, as.pop = "sex")
  sex <- glsex@other$ind.metrics$sex
  sexTidy <- factor(ifelse(grepl("f", sex, ignore.case = T), "f", "m"))
  glsex@other$ind.metrics$sex <- sexTidy
  
  glfilter <- glsex
  
  #number of populations
  n.pop= subpops
  #popSize = 20
  #number of generations to run
  n.gen= 20
  
  m.rate = 1e-7 #mutation rate
  if(dispersalType == "positive") mixev =  c(5,6,15,16)#mix generatoion
  if (dispersalType == "negative") mixev = c(1:20)[!(1:20 %in% c(5,6,15,16))]
  if (dispersalType == "constant") mixev = c(1:20)
  pop.inc = c(5,15,21) 
  m.ind = dispersal #mixing ratio - just going with 1 ind
  calcev = 1 #calc fst every calcev generation DONT NEED
  
  #start with a list of population based on pstart/glBase
  p <- list()
  pstart <- glfilter
  
  
  for (j in 1:n.pop){
    dads <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "m")
    mums <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "f")
    
    p[[j]] <- gl.sim.offspring(fathers = dads, mothers = mums,
                               noffpermother = 4, sexratio = 0.5)
    n.ind <- nInd(p[[j]])
    p[[j]]@pop <- factor(rep(paste0("P", j), n.ind))
    
    indNames(p[[j]]) <- paste0("p",j,"_",1:n.ind, "_",sprintf("%02d", 0))
  } #end for loop
  
  
  #generation zero
  em.gl.indmetrics <- function(glsim, gen = 0, disp.rate = 0){
    
    genn <- paste0("gen", sprintf("%02d", gen))
    glsim@other$ind.metrics <- data.frame(id = glsim@ind.names,
                                          pop = glsim@pop,
                                          trip = genn,
                                          dispersal = disp.rate, 
                                          sex = glsim@other$sex)
    return(glsim@other$ind.metrics)
    
  }
  
  glList <- list()
  pIndMet <- do.call("rbind", lapply(p, em.gl.indmetrics))
  pp <- do.call("rbind", p)
  pp@other$ind.metrics <- pIndMet
  
  cc <- 1
  glList[[cc]] <- pp 
  
  
  #number of mixing combinations between pops
  pcomb <- t(combn(n.pop,2))
  
  
  for(i in 1:n.gen){
    
    
    if ((i %in% mixev)){ #mix every mixev generations
      
      for (iii in 1:nrow(pcomb)){
        pa <- pcomb[iii,1]
        pb <- pcomb[iii,2]
        
        
        n.inda <- nInd(p[[pa]]) # different number of individuals
        n.indb <- nInd(p[[pb]]) # per population
        
        nmi <- ceiling(min(n.indb, n.inda)*m.ind) #number of individuals that swap
        
        pam <- sample(1:n.inda, nmi, replace = F)
        pbm <- sample(1:n.indb, nmi, replace = F)
        pback <- rbind(p[[pa]][-pam,], p[[pb]][pbm,])
        
        
        sexA <- factor(c(as.character(p[[pa]][-pam,]@other$sex),
                         as.character(p[[pb]][pbm,]@other$sex))) 
        pback@other$sex <- sexA                          
        popNames(pback) <- c(paste0("P",pa),paste0("P",pa))
        
        
        
        pback2 <- rbind(p[[pb]][-pbm,], p[[pa]][pam,]) 
        
        sexB <- factor(c(as.character(p[[pb]][-pbm,]@other$sex), 
                         as.character(p[[pa]][pam,]@other$sex)))
        pback2@other$sex <- sexB
        
        
        popNames(pback2) <- c(paste0("P",pb),paste0("P",pb))
        p[[pa]] <- pback
        p[[pb]] <- pback2
      }
    }
    
    
    for (ii in 1:n.pop){ # they breed creating the next generation
      
      gen <- p[[ii]]
      
      genDad <- gen[gen@other$sex == "male",]
      
      genMum <- gen[gen@other$sex == "female",]
      
      gen <- gl.sim.offspring(fathers = genDad, mothers = genMum,
                              noffpermother = nOff, sexratio = 0.5)
      
      genF <- gen[gen@other$sex == "female",]
      genF25 <- genF[sample(1:nInd(genF), min(nInd(genF), popSize/2)),]
      if ((i %in% pop.inc | (i-1) %in% pop.inc)) genF25 <- genF
      
      genM <- gen[gen@other$sex == "male",]
      genM25 <- genM[sample(1:nInd(genM), min(nInd(genM), popSize/2)),]
      if ((i %in% pop.inc | (i-1) %in% pop.inc)) genM25 <- genM
      
      
      gen50 <- rbind(genM25, genF25)
      gen50@other$sex <- factor(c(as.character(genM25$other$sex),
                                  as.character(genF25$other$sex)))
      
      
      n.ind <- nInd(gen50)
      gen50@pop <- factor(rep(paste0("P", ii), n.ind))
      indNames(gen50) <- paste0("p",ii,"_",
                                sprintf("%02d", 1:n.ind),
                                "_",sprintf("%02d", i))
      
      p[[ii]] <- gen50
      
      
    }
    
    # save generation
    
    
    dr <- ifelse(i %in% mixev, dispersal, 0)
    p2 <- lapply(p, function(x) em.gl.indmetrics(x, gen = i, disp.rate = dr))
    pIndMet <- do.call("rbind", p2)
    pp <- do.call("rbind", p)
    pp@other$ind.metrics <- pIndMet
    
    cc <- cc+1
    glList[[cc]] <- pp 
    
    
    
    cat(paste("Generation:", i,"/" ,n.gen, "\n"))
    flush.console()
    
    
  }
  
  
  indmet <- do.call("rbind", lapply(glList, function(x) x@other$ind.metrics))
  res <- do.call("rbind", glList)
  res@other$ind.metrics <- indmet
  res@other$loc.metrics.flags <- glBase@other$loc.metrics.flags
  
  
  return(res)
  
}


#' wrapper function for em.sim.boombust.ddd
#' 
#' This function simulates the three different types of dispersal on one
#' genlight object with a select dispersal rate and number of subpopulations.
#' Population size of subpopulations is set at 20. 
#' 
#' 
#' @param glBase -- a genlight object to base simulation on (starting snps)
#' @param subpops -- number of subpopulations
#' @param dispersal -- amount of dispersal
#' @return returns genlight for 20 generations
#' @export
#' @import tidyverse
#' @import dartR
#' @author Emily Stringer
#' @examples sim <- em.sim.boombust.ddd(gl, supops = 100, dispersal = 0.05)



em.sim.wrapper <- function(glBase, subpops = 25, dispersal = 0.10) {
  #positive
  simP <- em.sim.boombust.ddd(glBase, dispersal = 0.05,
                              dispersalType = "positive")
  simList <- list()
  disp <- c("positive", "negative", "constant")
  simList <- lapply(disp, function(x) em.sim.boombust.ddd(glBase, 
                                                          subpops,
                                                          dispersal,
                                                          dispersalType = x))
  names(simList) <- disp
  return(simList)
  
}
