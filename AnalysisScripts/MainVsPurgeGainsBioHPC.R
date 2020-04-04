library(BreedingSchemeLanguage)
if (exists("simEnv")){
  rm(list=names(simEnv), envir=simEnv)
  rm(simEnv)
}

# Set correct folders depending on compute environment
bioHPC <- TRUE
if (bioHPC){
  setwd("/home/jj332/PurgeDataScripts/")
  prefix <- "/home/jj332/PurgeDataScripts/"
} else{
  setwd("~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/PurgingSimulationScript2019")
  prefix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/"
}

allOutfiles <- list.files(paste(prefix, "AllOut", sep=""))
for (aof in allOutfiles){
  
  fileName <- paste(prefix, "AllOut/", aof, sep="")
  mainPurge <- readRDS(fileName)
  
  simEnv <- new.env()
  simEnv$sims <- mainPurge
  simEnv$nSim <- length(mainPurge)
  simEnv$nCore <- 1
  # Vectors of population means relevant for different selection events
  # The first population mean is always zero and consistent across all vectors
  # mainVals: main (=outcrossed) base population, main selected population
  mv <- plotData(popID = list(0:1, 0:1, 1, c(4:5, 13, 19), 5, c(8:9, 18), 9, c(12,14), 14, c(17,20), 20))
  # purgeVals: selfed base population, selfed selected population
  pv <- plotData(popID = list(0:1, 2:3, 3, 6:7, 7, 10:11, 11, 15:16, 16, 21:22, 22))
  # lossVals: selfed selected population, random-mated population from it
  lv <- plotData(popID = list(0:1, 3, c(4:5, 13, 19), 7, c(8:9, 18), 11, c(12,14), 16, c(17,20), 22, 23))
  # selfVals: outcrossed selected population, selfed population from it
  sv <- plotData(popID = list(0:1, 1, 2:3, 5, 6:7, 9, 10:11, 14, 15:16, 20, 21:22))
  saveRDS(list(mainVals=mv$data, purgeVals=pv$data, lossVals=lv$data, selfVals=sv$data), file=paste(prefix, "MainPurge/data", aof, sep=""))
  
  if (exists("simEnv")){
    rm(list=names(simEnv), envir=simEnv)
    rm(simEnv)
  }
  rm(list=setdiff(ls(), c("prefix", "forSomo", "allOutFiles"))); gc()
}
