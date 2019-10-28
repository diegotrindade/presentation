library(here)
mainDir<-here("Presentations")
setwd(mainDir)
source(here("Presentations", "GeneralPARAMETERS.R"))
source(here("Presentations", "SimulationFunctions_Suitab.R"))

# 
# subDir<-paste0("Different SAD Scenarios")
# dir.create(file.path(mainDir, subDir))
# setwd(file.path(mainDir, subDir))
# mainDir<-getwd()





#Load functions for the simulations (dispersal, competition...),
#and the parameters used (number of species and communities, etc.)


spRegion <- matrix(NA, nrow=nSpRegion, ncol=2,
                   dimnames = list(paste0("Sp.", 1:nSpRegion),
                                   c("Optima", "Disp")))
spRegion[ , "Optima"] <- rnorm(nSpRegion, 0, sdHabitat)
#spRegion[ , "Optima"] <- runif(nSpRegion, -sdHabitat, sdHabitat)

#spRegion[, "Disp"] <-runif(nrow(spRegion), 0,1) #In case you want to have different dispersal abilities
spRegion[, "Disp"] <-0.5 # This gives all sps the same dispersal ability

spHabitats<-spRegion[ , "Optima"]


### SITES AND ENVIRONMENTS
cells <- commID <- matrix(NA, ncol=sqrt(totalNCells), nrow=sqrt(totalNCells),
                          dimnames = list(paste0("R.",1:sqrt(totalNCells)),
                                          paste0("C.", 1:sqrt(totalNCells))))
ID<-1
#estava runif
HabitatsComm<-runif(sizeRegion^2, -sdHabitat, sdHabitat) #Species environment is drawn from a uniform distrib
#HabitatsComm<-rnorm(sizeRegion^2, 0, sdHabitat) #Species environment is drawn from a uniform distrib
for(i in 0:(sizeRegion-1)) {#Sites in rows
  for (j in 0:(sizeRegion-1)){
    #HabitatsAux<-rnorm(1, 0, sdHabitat)
    columsSelMax <- ((j+1)*sizeComm)
    cells[((i*sizeComm)+1):((i+1)*sizeComm) , (columsSelMax-sizeComm+1):columsSelMax] <- HabitatsComm[ID]
    commID[((i*sizeComm)+1):((i+1)*sizeComm) , (columsSelMax-sizeComm+1):columsSelMax] <- paste0("Comm.", formatC(ID, width = 3, format = "d", flag = "0"))
    ID<-ID+1
  }
}

#range(HabitatsCommChange)

par(mfrow = c(2,2))

# hist(HabitatsComm)
# hist(HabitatsComm+5)
# hist(spHabitats)
#####################
#####################
#####################SIMULATIONS
#####################
#####################



nTimeSteps<-200
envChangeTime<-100
encChangeSize<-5
HabitatsCommChange <- HabitatsComm + encChangeSize
saveSteps<-50
library(RColorBrewer)
library(lme4)
library(mgcv)
assymetryPar <- 10 #Parameter controlling competition (high assymetry gives more advantage to the species closest to the local environment)
mortalityPar <- 0.0001 #Probability that an individual dies in a given time step
immigrationPar <- 0.1 #Probability that an individual from the regional pool arrives in a given time step in each cell
reproductionPar <- 0.4
rangeEstablishment <- 3
suitThres <- 0.3
circular<-F #Shall the environment be "circular" (high and low extreme values will then be close to each other)
scenarios<-list()
scenarios[[1]] <- c(assymetry=assymetryPar, mortality=mortalityPar,  
                    immigration=immigrationPar, establishment=rangeEstablishment,
                    reproduction = reproductionPar)
names(scenarios)[1] <- paste0("Random." , 6)


for(s in 1:length(scenarios)){
  cat(paste("\n\n\n\n Starting with scenario =", names(scenarios)[s], "\n\n\n\n"))
  subDir<-"Uniform Distrib_Env"
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  commMat <- matrix(NA, ncol=ncol(cells), nrow=nrow(cells),
                    dimnames=dimnames(cells))
  par(mfrow=c(1,2), mar=c(1,1,1,1))
  SaveID<-1
  startTime<-1
  
  
  observedRichness <- matrix(NA, nrow=length(unique(as.vector(commID))), ncol= nTimeSteps,
                             dimnames = list(paste0("Comm.", 1:length(unique(as.vector(commID)))),
                                             paste0("Time.", 1:nTimeSteps)))
  establishRichness <- matrix(NA, nrow=length(unique(as.vector(commID))), ncol= nTimeSteps,
                              dimnames = list(paste0("Comm.", 1:length(unique(as.vector(commID)))),
                                              paste0("Time.", 1:nTimeSteps)))
  exctDebtRichness <- matrix(NA, nrow=length(unique(as.vector(commID))), ncol= nTimeSteps,
                             dimnames = list(paste0("Comm.", 1:length(unique(as.vector(commID)))),
                                             paste0("Time.", 1:nTimeSteps)))
  poolSize <- matrix(NA, nrow=length(unique(as.vector(commID))), ncol= nTimeSteps,
                     dimnames = list(paste0("Comm.", 1:length(unique(as.vector(commID)))),
                                     paste0("Time.", 1:nTimeSteps)))
  darkDiv <- matrix(NA, nrow=length(unique(as.vector(commID))), ncol= nTimeSteps,
                    dimnames = list(paste0("Comm.", 1:length(unique(as.vector(commID)))),
                                    paste0("Time.", 1:nTimeSteps)))
  for (time in startTime:nTimeSteps){
    #########PARAMETERS THAT CHANGE ALONG WITH ENVIRONMENT
    if(time == envChangeTime){ #A sudden environmental change (increase everywhere)
      cells <- cells + encChangeSize
      HabitatsCommChange <- HabitatsComm + encChangeSize
    }
    
    ####which species can establish in each community when time==1?
    # establishBefore <- establishAfter <-  matrix(NA, nrow=nrow(spRegion), ncol=length(HabitatsComm),
    #                           dimnames = list(rownames(spRegion), paste0("Comm.", 1:length(HabitatsComm))))
    if (time == 1){
      suitabilityMat <- suitability (cells=cells, spHabitats=spHabitats, 
                                     sd_range = scenarios[[s]]["establishment"],
                                     circular=circular)
      #Lets give suitability = 0 to species with suitability < 0.05 (or another vakue)
      suitabilityMat <- replace(suitabilityMat, suitabilityMat < suitThres, 0)
    }
    # HabitatsComm=HabitatsComm+3
    if (time == (envChangeTime + 1)){
      suitabilityMat <- suitability (cells=cells, spHabitats=spHabitats, 
                                     sd_range = scenarios[[s]]["establishment"],
                                     circular=circular)
      #Lets give suitability = 0 to species with suitability < 0.05 (or another vakue)
      suitabilityMat <- replace(suitabilityMat, suitabilityMat < suitThres, 0)
      
    }
    
    ##############################################################
    if(time==1 | time%%saveSteps==0){
      if(time>1) {cat(paste("Last rep took", round(as.numeric(endTime-startTime),2), "secs\n"))}
      cat(paste("Step", time, "out of", nTimeSteps, "\n"))
    }
    startTime <- Sys.time()
    dispersionResult <- dispersal(cells=cells, localSp=commMat, spRegion=spRegion,
                                  probImmig = scenarios[[s]]["immigration"] , 
                                  probDisp = scenarios[[s]]["reproduction"])
    
    establishmentResult <- establishment(cells=cells, arrivingSp=dispersionResult, 
                                         spHabitats=spHabitats, suitability=suitabilityMat)
    
    competitionResult <- competition (localSp=commMat,
                                      establishedSp=establishmentResult,
                                      suitability=suitabilityMat, advLocal = 1,
                                      assymetry = scenarios[[s]]["assymetry"])
    
    if(time==1){
      
      dens<-numeric(length(cells))
    } else{
      dens<-numeric()
      for(i in 1:length(commID)){
        dens[i]<- comMatSim[commID[i], commMat[i]]
        dens[is.na(dens)] <- 0
      }
    }
    
    mortalityResult <- mortalityDensityDependent(cells=cells, localSp=competitionResult, density=dens,
                                                 probMortality = scenarios[[s]]["mortality"], 
                                                 maxMortalityInc = 0.3)
    
    
    mortalityResult <- mortalityDensityDependent(cells=cells, localSp=competitionResult, density=dens,
                                                 probMortality = scenarios[[s]]["mortality"], 
                                                 maxMortalityInc = 0.3)
    commMat <- mortalityResult
    
    
    #Plot how it is going from time to time (each 100 steps)
    if(time%%10==0 |time==1){
      x<-1:nrow(cells)
      y<-1:ncol(cells)
      image(x, y, matrix(spHabitats[commMat], ncol=max(y), nrow=max(x)), col = brewer.pal(11, "Spectral"), axes = FALSE )
      abline(h=seq(0,max(x), by=sizeComm))
      abline(v=seq(0,max(y), by=sizeComm))
      image(x, y, matrix(cells, ncol=max(y), nrow=max(x)), col = brewer.pal(11, "Spectral"), axes = FALSE )
      abline(h=seq(0,max(x), by=sizeComm))
      abline(v=seq(0,max(y), by=sizeComm))
    }
    #Save results when it is time to save
    
    # if(time>burnIn & (time-burnIn)%%saveSteps==0){ ###SAVING
    #   cat(paste("Saving results after", time, "reps\n"))
    commIDVect<-as.vector(commID)
    spComp<-tapply(commMat, commIDVect, table)
    
    envOptSpComp<-tapply(spRegion[,1][commMat], commIDVect, unique, na.rm=T)
    comMatSim<-matrix(0, nrow=length(unique(as.vector(commID))), ncol= nSpRegion,
                      dimnames = list(names(spComp), rownames(spRegion)))
    for(j in 1:length(spComp)){
      if(length(spComp[[j]]>0)){
        spSelAux<-paste0("Sp.",names(spComp[[j]]))
        comMatSim[j, spSelAux] <- spComp[[j]]
      }
    }
    
    observedRichness[ , time] <- rowSums(comMatSim>0)
    
    obsAux<-establishAux <- exctDebtRichnessAux <- poolSizeAux <- darkDivAux <- numeric(length(HabitatsCommChange))
    
    minEstablish <- HabitatsCommChange - rangeEstablishment
    maxEstablish <- HabitatsCommChange + rangeEstablishment
    
    for( comm in 1:length(HabitatsCommChange)){
      IDStablished <- which(suitabilityMat[comm, ] > 0)
      establishAux[comm] <- length(IDStablished)
      exctDebtRichnessAux[comm] <- sum(comMatSim[comm, -IDStablished]>0)
      obsAux[comm] <- sum(comMatSim[comm, ])
      poolSizeAux[comm] <- length(IDStablished)
      darkDivAux[comm] <- poolSizeAux[comm] - obsAux[comm]
      
      #darkDivAux[comm] <- sum(comMatSim[comm, IDStablished]>0)
      
    }
    
    establishRichness[, time] <- establishAux
    exctDebtRichness[, time] <- exctDebtRichnessAux
    darkDiv[, time] <- darkDivAux
    poolSize[, time] <- poolSizeAux
    
    darkDiv[,time]<-  poolSize[, time] - observedRichness[ , time]
    
    poolSize[, time] <- observedRichness[ , time] + darkDiv[, time]
    
    
    # write.table(comMatSim, file=paste0("ResultsSim.", time,".txt"), quote=F)
    # }
    
    
    
    endTime <- Sys.time()
  }
}



jpeg(here("Presentations", file =  "ind_basee.png"), width = 7, height = 5, units = 'in', res = 300)

#png(here("Presentations", file =  "ind_basee.png", res=144))
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:25){
  plot(1:time, observedRichness[i,], type="l", main=sort(HabitatsCommChange, decreasing = T)[i], col = "blue")
  abline(v=100, lty=2, col=2)
  par(new=TRUE)
  plot(1:time, darkDiv[i,], type="l", yaxt='n', ann=FALSE, lth=3)
  #par(new=T)
  #plot(1:time, exctDebtRichness[i,], type="l", col="red")
}
#library(ggplot2)

dev.off()


?ggsave

for(i in 1:25){
  plot(1:time, darkDiv[i,], type="l", main=sort(HabitatsCommChange, decreasing = T)[i])
  abline(v=250, lty=2, col=2)
}

for(i in 1:25){
  plot(1:time, poolSize[i,], type="l", main=sort(HabitatsComm, decreasing = T)[i])
  abline(v=250, lty=2, col=2)
}

jpeg(here("Presentations", file =  "extd.png"), width = 7, height = 5, units = 'in', res = 300)

par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:25){
  plot(1:time, exctDebtRichness[i,],col="red", type="l", main=sort(HabitatsCommChange, decreasing = T)[i])
  abline(v=250, lty=2, col=2)
  
}


dev.off()

range(HabitatsCommChange)

par(mfrow=c(2,2), mar=c(2,2,2,2))
hist(HabitatsCommChange)
hist(spHabitats)

max(HabitatsComm)
mean(spHabitats)
range(HabitatsCommChange)
sd(HabitatsComm-5)
range(spHabitats)