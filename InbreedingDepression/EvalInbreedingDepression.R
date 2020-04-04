# Set correct folders depending on compute environment
library(BreedingSchemeLanguage)

bioHPC <- FALSE
if (bioHPC){
  setwd("/home/jj332/PurgeDataScripts/")
  prefix <- "/home/jj332/PurgeDataScripts/"
} else{
  setwd("~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/PurgingSimulationScript2019")
  prefix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/"
}

allOutFiles <- list.files(path=paste0(prefix, "AllOut"), pattern="_purge")

print("####### outFiles #######")
print(paste0(prefix, "AllOut"))
print(allOutFiles)

for (fileName in allOutFiles){
  print(fileName)
  allOut <- readRDS(paste0(prefix, "AllOut/", fileName))
  simEnv <- new.env()
  simEnv$sims <- allOut
  simEnv$nSim <- length(allOut)
  simEnv$nCore <- 1
  simEnv$onlyCost <- FALSE

  outSuff <- substring(fileName, 30)
  
  selfFertilize(nProgeny = 200, popID = 0:1) # create population 24
  selfFertilize(nProgeny = 200, popID = 23)  # create population 25
  inbreedStrtEnd <- plotData(popID=list(0:1, 24, 23, 25))
  saveRDS(inbreedStrtEnd, file=paste0(prefix, "Inbred/inbredData_", outSuff))
  
  inbreedStrtEndcsv <- as.data.frame(inbreedStrtEnd$data)
  write.csv(inbreedStrtEndcsv, file=paste0(prefix, "Inbred/inbredCSV_", outSuff, ".csv"))

  if (exists("simEnv")){
    rm(list=names(simEnv), envir=simEnv)
    rm(simEnv)
  }
  rm(list=setdiff(ls(), c("allOutFiles", "fileName", "prefix"))); gc()
}

setwd("./InbreedingDepression")
inbrFiles <- list.files(path=".", pattern=".csv")
inbrFiles <- inbrFiles[-grep("Additive", inbrFiles)]
print("####### inbrFiles #######")
print(inbrFiles)

allInbrDep <- data.frame()
for (fileName in inbrFiles){
  inbrDF <- read.csv(fileName)
  inbrDF <- inbrDF[inbrDF$group != 25,]
  foundSize <- as.numeric(strsplit(fileName, split="[FN]")[[1]][2])
  scheme <- strsplit(fileName, split="[_A]")[[1]][4]
  cycle <- rep(c(0, 5), 24)
  replication <- rep(1:24, each=2)
  inbrDep <- diff(inbrDF$g)[1:48 * 2 - 1]
  allInbrDep <- rbind(allInbrDep, data.frame(scheme=scheme, foundSize=foundSize, cycle=cycle, replication=as.factor(replication), inbrDep=-inbrDep))
}

fitFnd <- lm(inbrDep ~ replication + foundSize*scheme, data=allInbrDep[allInbrDep$cycle == 0,])
summary(fitFnd)
anova(fitFnd)

fitC5 <- lm(inbrDep ~ replication + foundSize*scheme, data=allInbrDep[allInbrDep$cycle == 5,])
summary(fitC5)
anova(fitC5)

allInbrDep$sizeSchem <- paste0(allInbrDep$scheme, allInbrDep$foundSize)

fitFnd <- lm(inbrDep ~ replication + sizeSchem, data=allInbrDep[allInbrDep$cycle == 0,])
summary(fitFnd)
anova(fitFnd)
c0means <- coefficients(fitFnd)
names(c0means)[1] <- "sizeSchemMeanSelf1500"
c0grndMn <- c0means[1] + sum(c0means[grep("repl", names(c0means))])/24
c0means <- c0means[-grep("repl", names(c0means))]
c0means[2:9] <- c0means[2:9] + c0grndMn

fitC5 <- lm(inbrDep ~ replication + sizeSchem, data=allInbrDep[allInbrDep$cycle == 5,])
summary(fitC5)
anova(fitC5)
c5means <- coefficients(fitC5)
names(c5means)[1] <- "sizeSchemMeanSelf1500"
c5grndMn <- c5means[1] + sum(c5means[grep("repl", names(c5means))])/24
c5means <- c5means[-grep("repl", names(c5means))]
c5means[2:9] <- c5means[2:9] + c5means[1]
