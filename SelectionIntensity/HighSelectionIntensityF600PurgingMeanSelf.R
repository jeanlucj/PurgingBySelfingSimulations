#Repeat Simulation in April, 2019 for HighSelectionIntensity=SelfMean
library(rrBLUP)
library(snowfall)
library(BreedingSchemeLanguage)
# Mean phenotype based on maternal parent
# Makes simplifying assumption of uniform error variances
makeSynPheno <- function(sEnv=simEnv, pID){
  parent.env(sEnv) <- environment()
  makeSynPheno.func <- function(bsl, pID){
    gr <- subset(bsl$genoRec, popID %in% pID)
    pr <- subset(bsl$phenoRec, phenoGID %in% gr$GID)
    seedPar <- sort(unique(gr$pedigree.1))
    makeSyn <- function(sp){
      progGID <- gr$GID[gr$pedigree.1 %in% sp]
      spr <- subset(pr, phenoGID %in% progGID)
      res <- colMeans(spr[, c("error", "pValue")])
      res[1] <- res[1] / nrow(spr)
      return(res)
    }
    synPhen <- t(sapply(seedPar, makeSyn))
    toAdd <- data.frame(phenoGID=seedPar, loc=1, year=1, synPhen)
    colnames(toAdd) <- c("phenoGID", "loc", "year", "error", "pValue")
    bsl$phenoRec <- rbind(bsl$phenoRec, toAdd)
    return(bsl)
  }
  with(sEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, makeSynPheno.func, pID=pID)
      sfStop()
    }else{
      sims <- lapply(sims, makeSynPheno.func, pID=pID)
    }
  })
}

# How many generations of evolution
nGen <- "2000"
# Decide how many haplotypes, markers, and qtl
oracle <- FALSE
selMrk <- FALSE
slimPopSize <- 500
nHap <- 2*slimPopSize
nMrk <- 1200
nQTL <- 300
# Set up matrix showing how many simulations and how many chromosomes per sim
nSim <- 24
nChr <- 5

baseNameDef <- paste("_purgeG", nGen, "H", nHap, "M", nMrk, "Q", nQTL, "N", slimPopSize, ifelse(oracle, "ORA", "NOR"), ifelse(selMrk, "SEL", "RAN"), sep="")

load(file=paste(baseNameDef, "DefSpec.RData", sep=""))
simEnv <- new.env()
redSim <- 24
#simEnv$sims <- sims
simEnv$sims <- sims[1:redSim]#for shorter time
simEnv$nSim <- nSim
simEnv$nCore <- 2
simEnv$onlyCost <- FALSE
rm(sims); gc()

#HighSelectionIntensity=30=5%
nFound <- 600
nSelect <- 30
nProgeny <- 1000
toCET <- nProgeny / 2 - nSelect
toPYT <- nProgeny / 4 - nSelect
finalPop <- 200
nSelfed <- 10
nSelWF <- 5

baseName <- paste("_purgeG", nGen, "H", nHap, "M", nMrk, "Q", nQTL, "N", slimPopSize, "F", nFound, ifelse(oracle, "ORA", "NOR"), ifelse(selMrk, "SEL", "RAN"), sep="")

# Start with more cores because the problem is not too big yet
s <- Sys.time()
###############################################################
library(BreedingSchemeLanguage)
defineVariances(gByLocVar=0, gByYearVar=0, plotTypeErrVars=c(Founder=1, PYT=4, CET=9, Seedling=16))

# Phenotypic selection here
initializePopulation(nInd=nFound) # popID 0
genotype(popID=0)
phenotype(plotType="Founder")
predictValue(popID=0, sharingInfo="markers")
select(nSelect=nSelect, popID=0) # popID 1

selfFertilize(nProgeny=nSelect*nSelfed, popID=1) # popID 2
phenotype(plotType="Seedling", popID=2)
select(nSelect=nSelWF, random=F, type="WithinFamily", popID=2) # popID 3
makeSynPheno(pID=2:3)
cross(nProgeny=nProgeny, notWithinFam=T) # popID 4
genotype(popID=4)

print(paste("Cycle 1 Created", Sys.time() - s))
predictValue(popID=4, sharingInfo="markers")
select(nSelect=nSelect, popID=4) # popID 5

selfFertilize(nProgeny=nSelect*nSelfed, popID=5) # popID 6
phenotype(plotType="Seedling", popID=6)
select(nSelect=nSelWF, random=F, type="WithinFamily", popID=6) # popID 7
makeSynPheno(pID=6:7)
cross(nProgeny=nProgeny, notWithinFam=T) # popID 8
genotype(popID=8)

print(paste("Cycle 2 Created", Sys.time() - s))
predictValue(popID=8, sharingInfo="markers")
select(nSelect=nSelect, popID=8) # popID 9

selfFertilize(nProgeny=nSelect*nSelfed, popID=9) # popID 10
phenotype(plotType="Seedling", popID=10)
select(nSelect=nSelWF, random=F, type="WithinFamily", popID=10) # popID 11
makeSynPheno(pID=10:11)
cross(nProgeny=nProgeny, notWithinFam=T) # popID 12
genotype(popID=12)

# Reduce number of cores to be on the safe side
simEnv$nCore <- 1

print(paste("Cycle 3 Created", Sys.time() - s))
select(nSelect=toCET, popID=4, random=T) # popID 13
phenotype(popID=c(5, 13), plotType="CET")

predictValue(popID=12, sharingInfo="markers")
select(nSelect=nSelect, popID=12) # popID 14

selfFertilize(nProgeny=nSelect*nSelfed, popID=14) # popID 15
phenotype(plotType="Seedling", popID=15)
select(nSelect=nSelWF, random=F, type="WithinFamily", popID=15) # popID 16
makeSynPheno(pID=15:16)
cross(nProgeny=nProgeny, notWithinFam=T) # popID 17
genotype(popID=17)

print(paste("Cycle 4 Created", Sys.time() - s))
select(nSelect=toCET, popID=8, random=T) # popID 18
phenotype(popID=c(9, 18), plotType="CET")
select(nSelect=toPYT, popID=13, random=T) # popID 19
phenotype(popID=c(5, 19), plotType="PYT")

predictValue(popID=17, sharingInfo="markers")
select(nSelect=nSelect, popID=17) # popID 20

selfFertilize(nProgeny=nSelect*nSelfed, popID=20) # popID 21
phenotype(plotType="Seedling", popID=21)
select(nSelect=nSelWF, random=F, type="WithinFamily", popID=21) # popID 22
makeSynPheno(pID=21:22)
cross(nProgeny=finalPop, notWithinFam=T) # popID 23
print(paste("Cycle 5 Created", Sys.time() - s))

outputResults(summarize=FALSE, saveDataFileName=paste(baseName, "HighSelectionIntensityF600_MeanSelfAllOut", sep=""))
pdf(paste(baseName, "HighSelectionIntensityF600_MeanSelfPlot.pdf", sep=""))
plotData(addDataFileName=paste(baseName, "HighSelectionIntensityF600_MeanSelfPlotOut", sep=""), popID=list(0:1, c(4:5,13,19), c(8:9,18), c(12,14), c(17,20), 23))
dev.off()
