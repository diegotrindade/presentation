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

set.seed(2)
nSpRegion <- 100
sizeRegion <- 20 # the total number of communities is sizeRegion*sizeRegion
sizeComm <- 5 # the total number of cells in each comm is sizeComm*sizeComm
totalNCells <- sizeRegion^2 * sizeComm^2
#lengthGenome <- 6#Number of genes considered
nHabitats <- 12 #Number of different environments
Habitats <- 1:nHabitats
freqHabitats <- "equal" # Either equal or a relative freq for each habitat
probSpInHabitats <- "equal" # Either equal or a relative freq for each habitat within species genomes
competition <- "rel" # c( "rel", "abs"). "rel" means competition is estimated as a probability/ "abs" means the stronger takes it all
localAdv <- 1.5 # How much does the relative comp ability should be increased? 
sdHabitat<-10

