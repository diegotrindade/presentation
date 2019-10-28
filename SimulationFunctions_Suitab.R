### Metacommunity formed by a series of sites (30*30?), each one formed by a
### series of cells (10x10). All cells start being empty at the beggining, and 
### we will simulate community dynamics including four driving processes (Vellend):
### A) Dispersal: In each cell, species will be able to arrive randomly from a
###     regional pool. From this mechanism, all species have equal probability of
###     arrival to any site. In addition, species will also be able to disperse
###     to any cell. Established individuals will be randomly chosen to disperse
###     or not in each time step (lets say 10% of the time). Species will have 
###     different (or not) dispersal kernels, which will be given by some genetic
###     characteristic (related or not to other abilities).
### B) Drift: There is a mortality rate, which can be habitat specific (or more
###     complex. such as habitat*species, or whatever). In each tme step, 
###     the individual living in a cell will have a probability to die that is 
###     given by this rate.
### C) Selection: Individuals will be more or less adapted to specific habitats
###     based on their genomes. Each individual will have a genome composed of a
###     number of letters (genes?). Each letter correspond to the ability to live
###     in a given environment; the higher the proportions of those letters in 
###     the genome, the higher the competitive ability of that species in that
###     environment. When a species disperse into a cell, if there is another
###     one living there, they will compete. The winner will stay in the cell.
###     (How to do this? I think it makes more sense as a ratio
###     between the competitive abilities of the species, perhaps giving more of 
###     an advantage to the species that is occupying the cell already. This ratio
###     will give us a probability for each species, and we sample from that. The
###     alternative: the stronger wins always, regardless of the difference).
### D) Speciation: not used here (although it kind of determines competitive
###     and dispersal abilities)   

######
######
### 0) SUITABILITY
# Suitability indicates how cloes an environment is to the optimum of a certain
# species. In general, if the environment value is equal to the species, then 
# suitability will be 1, and decrease towards 0 as distance increases (following
# a normal distribution). sd_range gives us the sd of the normal distribution.
suitability <- function (cells, spHabitats, sd_range = 1, circular=F){
  commIDVect<-as.vector(commID)
  EnvCommuVect<-as.vector(cells)
  EnvComm<-tapply(EnvCommuVect, commIDVect, mean)
  suitabilityMatrix <- matrix(NA, nrow=length(EnvComm), ncol=length(spHabitats),
                              dimnames = list(names(EnvComm), names(spHabitats)))
  rangeEnv <- range(range(cells), range(spHabitats))
  lengthEnv<- abs(rangeEnv[2] - rangeEnv[1])
  
  for(i in 1:nrow(suitabilityMatrix)){
    distancesAux <- abs(spHabitats - EnvComm[i])
    if(circular){
      distancesAux <- ifelse(distancesAux < lengthEnv/2, distancesAux, lengthEnv - distancesAux)
    }
    suitabilityMatrix[i,] <- distancesAux
  }
  #### Now we can transform suitability to probability according to a normal
  #### distribution based on the distance.
  suitabilityMatrix <- dnorm(suitabilityMatrix, mean=0, sd=sd_range) / 
    dnorm(0, mean=0, sd=sd_range)
  
  return(suitabilityMatrix)
} 



### SIMULATIONS
### 1) IMMIGRATION FROM LOCAL POOL
### 2) DISPERSAL (FOR NEXT TIME STEP)
# Immigration is a random process, with some prevalence
# dispersal is estimated as a lognormal distribution. Lets use the maximum length
# of our matrix as the reference for maximum dispersal distance:
# 
# 
dispersal <- function(cells, localSp, spRegion, probImmig = 0.05 , probDisp = 0.4){
  rowsCells<-nrow(localSp)
  colsCells<-ncol(localSp)
  #A) Immigration from the pool  
  immigrationEvent <- rbinom(n=length(cells), size=1, prob = probImmig)
  immigrationEvent <- replace(immigrationEvent, immigrationEvent==0, NA)
  immigrationEvent[which(immigrationEvent==1)] <- sample(1:nrow(spRegion), size=sum(immigrationEvent, na.rm=T), replace=T)
  immigrationEvent <- matrix(immigrationEvent, ncol=colsCells, nrow=rowsCells,
                             dimnames=dimnames(cells))
  #B) Dispersal from established plants  
  dispersalSP <- matrix(NA, ncol=colsCells, nrow=rowsCells,
                        dimnames=dimnames(cells))
  if(any(!is.na(localSp) & probDisp>0)){
    maxdist <- round(as.numeric(dist(rbind(c(1,1), c(rowsCells,colsCells)))))
    occupiedCells <- which(!is.na(localSp))
    
    dispersingSps <- localSp[occupiedCells]
    dispersionEvent <- rbinom(n=length(dispersingSps), size=1, prob = probDisp)
    dispDirAux <-runif(length(dispersingSps), 0, 2*pi)
    dispDistAux <- rlnorm(n=length(dispersingSps),meanlog=log(maxdist)*0.2*spRegion[dispersingSps, "Disp"])
    dispDistAux <- pmax(1, pmin(dispDistAux, maxdist))
    dispDirAux[which(dispDistAux==0)] <- NA
    dispDirAux[which(dispersionEvent==0)] <- NA
    dispersingSps[which(dispersionEvent==0)] <- NA
    changeX<- round(dispDistAux * cos(dispDirAux))
    changeY<- round(dispDistAux * sin(dispDirAux))
    dispersalEventOriginX <- ceiling(occupiedCells/rowsCells)
    dispersalEventOriginY <- occupiedCells%%rowsCells +1
    dispersalEventFinalY <- 1 + (dispersalEventOriginY + changeY)%%rowsCells
    dispersalEventFinalX <- 1 + (dispersalEventOriginX + changeX)%%colsCells
    dispersalEventFinalX[which(dispersionEvent==0)] <- NA
    dispersalEventFinalY[which(dispersionEvent==0)] <- NA
    dispersingSps[which(dispersionEvent==0)] <- NA
    for(i in 1:length(dispersingSps)){
      dispersalSP[dispersalEventFinalX[i], dispersalEventFinalY[i]]<-dispersingSps[i]
    }
  }
  #C) Final result (if two sps arrive to a place, one of them is chosen randomly)  
  dispersion <- pmax(immigrationEvent, dispersalSP, na.rm=T)
  coincidentCells <- which(!is.na(dispersalSP) & !is.na(immigrationEvent)) 
  randomWinner <- sample(c("IMM","DISP"), size=length(coincidentCells), replace=T)
  winner<-ifelse(randomWinner=="IMM", immigrationEvent[coincidentCells], dispersalSP[coincidentCells])
  dispersion[coincidentCells] <- winner
  return(dispersion)
}


### 3) Establishment
# Establishment simply evaluates if the species that arrive by any means are able to
# survive in a habitat for which they have low adaptation value (if a species optimum
# is within a cell's range, then it can establish there):
establishment <- function (cells, arrivingSp, spHabitats, suitability){
  establishmentSP <- arrivingSp
  arrivalCells <- which(!is.na(establishmentSP)) 
  arrivingSpID <- establishmentSP[arrivalCells] 
  commArriving <- commID[arrivalCells]
  commArrivingRows <- match(commArriving, rownames(suitability))
  probEstablishment <- suitability[cbind(commArrivingRows, arrivingSpID)]
  establishmentSPVector <- rbinom(n = length(probEstablishment), size=1, 
                                  prob = probEstablishment)
  establishmentSPVector <-replace(arrivingSpID, establishmentSPVector == 0, NA)
  establishmentSP[arrivalCells] <- establishmentSPVector
  
  return(establishmentSP) 
}

### 3) COMPETITION
### Competition plays a role when a propagule of a species arrives to an occupied cell.
### Each species has a degree of adaptation to a given habitat (0-1) given by its genome
### The probability that each species wins depends on this adaptation (given by the 
### ratio of suitabilities of the two species). Further, we can give advantage to 
### the species that was already there. 
competition <- function (localSp, establishedSp, suitability, advLocal = 1, 
                         assymetry = 1){
  winnerSp<-pmax(localSp, establishedSp, na.rm=T)
  coincidentCells <- which(!is.na(localSp) & !is.na(establishedSp)) 
  localSpVector <- localSp[coincidentCells]
  establishedSpVector <- establishedSp[coincidentCells]
  commCompetingVector <- commID[coincidentCells]
  commCompetingRows <- match(commCompetingVector, rownames(suitability))
  suitabilityLocal <- suitability[cbind(commCompetingRows, localSpVector)]
  suitabilityEstablished <- suitability[cbind(commCompetingRows, establishedSpVector)]
  compLocal <- suitabilityLocal * advLocal
  compEstablish <- suitabilityEstablished
  diffComp <- (compLocal - compEstablish)*assymetry
  probsLocalWin <- 1 / (1 + exp(-diffComp))
  localWinner <- rbinom(length(probsLocalWin), size=1, prob = probsLocalWin)
  winnerSpVector <- rep(NA, length(localWinner))
  winnerSpVector[which(localWinner==1)] <- localSpVector[which(localWinner==1)]
  winnerSpVector[which(localWinner==0)] <- establishedSpVector[which(localWinner==0)]
  winnerSp[coincidentCells] <- winnerSpVector
  
  return(winnerSp)
} 




### 4) MORTALITY
### Mortaility is simply a random process, which could be habitat-dependent if we
### want to. Each individual has a probability of dying.
### probMortality is either a global single number, or a probability for
### each habitat
mortalityDensityDependent <- function (cells, localSp, density, probMortality = 0.1, maxMortalityInc = 0.1){
  # max increase indicates how many times (in proportion) mortality increases compared to the baseline
  # for a species with max density
  maxDens <- sizeComm^2
  probMortalityCell <- probMortality + maxMortalityInc * dens/maxDens
  probMortalityCell[is.na(probMortalityCell)] <- 0
  death <- rbinom(n=length(cells), size=1, prob = probMortalityCell)
  localSp[which(death==1)] <- NA
  return(localSp)
} 


mortality <- function (cells, localSp, probMortality = 0.1){
  if(length(probMortality==1)){
    death <- rbinom(n=length(cells), size=1, prob = probMortality)
  } else{
    if(length(probMortality!=nHabitat)){ stop("Mortality must be either a single number or a vector with a mortality rate for each habitat")
    }
    mortalityVector<-probMortality[cells]
    death <- rbinom(n=length(cells), size=1, prob = mortalityVector)
  } 
  localSp[which(death==1)] <- NA
  return(localSp)
}