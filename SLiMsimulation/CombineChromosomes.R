# SLiM has made single chromosome simulations
# Concatenate multiple chromsomes into one genome
# Create a BSL sim environment with all of them to parallelize
# Run a BSL simulation

source("~/Documents/simOut/SLiM/parseSLiMfunctions.R")

# Make multi-chromosome environments
concatChr <- function(whichChr, hapVec, nMrk, nQTL, oracle=FALSE){
  mutData <<- readRDS(file=paste(basePath, "MutData", whichChr[1], ".rds", sep=""))
  indRaw <<- readRDS(file=paste(basePath, "IndRaw", whichChr[1], ".rds", sep=""))
  simEnv <- makeSimEnv(list(hapVec), nMrk, nQTL, oracle=oracle)
  cat("Chr started", whichChr[1], "\n")
  for (e in 2:length(whichChr)){
    mutData <<- readRDS(file=paste(basePath, "MutData", whichChr[e], ".rds", sep=""))
    indRaw <<- readRDS(file=paste(basePath, "IndRaw", whichChr[e], ".rds", sep=""))
    addEnv <- makeSimEnv(list(hapVec), nMrk, nQTL, oracle=oracle)
    nChr <- addSimEnvChr(simEnv, addEnv)
    cat("Chr concatenated", whichChr[e], "\n")
  }
  return(simEnv)
}

# How many generations of evolution
nGen <- "2000"
basePath <- paste("~/Documents/simOut/SLiM/consElem6Ne100out/constrainedElements6_", nGen, "Ne100Out", sep="")
slimProg <- readLines(paste("~/Documents/simOut/SLiM/constrainedElements6_", nGen, "HalfNeutMoreDomNe100.txt", sep=""))

# Decide how many haplotypes, markers, and qtl
oracle <- FALSE
slimPopSize <- 100
nHap <- 1000
nMrk <- 480
nQTL <- 120
# Set up matrix showing how many simulations and how many chromosomes per sim
nSim <- 24
nChr <- 5

defSpecName <- paste("_purgeG", nGen, "H", nHap, "M", nMrk, "Q", nQTL, "N", slimPopSize, ifelse(oracle, "ORA", "NOR"), "DefSpec.RData", sep="")

didMutData <- TRUE
if (!didMutData){
  # make mutData and indRaw objects for all of the simulated chromosomes
  for (whichSim in 0:17){
    slimOut <- readLines(paste(basePath, whichSim, ".txt", sep=""))
    mutData <- parseMut()
    saveRDS(mutData, file=paste(basePath, "MutData", whichSim, ".rds", sep=""))
    indRaw <- parseInd()
    saveRDS(indRaw, file=paste(basePath, "IndRaw", whichSim, ".rds", sep=""))
    cat("Made mutData and indRaw", whichSim, "\n")
  }
}

didCombineChr <- TRUE
if (!didCombineChr){
  chrToGenome <- sapply(1:nSim, function(d) sample(0:17, nChr))
  
  simEnv <- concatChr(chrToGenome[,1], hapVec=1:nHap, nMrk=nMrk, nQTL=nQTL, oracle=oracle)
  for (sim in 2:nSim){
    print("nSim combined so far")
    print(combineSimEnvSims(simEnv, concatChr(chrToGenome[,sim], hapVec=1:nHap, nMrk=nMrk, nQTL=nQTL, oracle=oracle)))
  }
  
  simEnv$nCore <- 6
 # _purgeG2000H1000M480Q120N100NORDefSpec.RData 
  with(simEnv, {
    save(sims, nSim, nCore, file=defSpecName)
  })
  
}

# This code should be used prior to the BSL simulation program.
if (FALSE){
  #else {
  load(file=defSpecName)
  simEnv <- new.env()
  simEnv$sims <- sims
  simEnv$nSim <- nSim
  simEnv$nCore <- nCore
  rm(sims)
  #}
}
