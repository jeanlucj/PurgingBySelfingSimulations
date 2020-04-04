# Set correct folders depending on compute environment
bioHPC <- FALSE
if (bioHPC){
  setwd("/home/jj332/PurgeDataScripts")
  prefix <- "/home/jj332/PurgeDataScripts/"
} else{
  setwd("~/Documents/BioHPC/jj332/PurgeDataScripts")
  prefix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/"
}

allOutFiles <- list.files(paste(prefix, "AllOut", sep=""))
# allOutFiles <- allOutFiles[-grep("PYT", allOutFiles)]
for (aof in allOutFiles){
  
  fileName <- paste(prefix, "AllOut/", aof, sep="")
  sims <- readRDS(fileName)
  
  popIDlist <- list(0:1, c(4:5, 13, 19), c(8:9, 18), c(12,14), c(17,20), 23)
  calcVars <- function(bslData){
    calcVar <- function(popID) return(var(bslData$gValue[bslData$genoRec$popID %in% popID]))
    return(sapply(popIDlist, calcVar))
  }
  varsByRepCyc <- sapply(sims, calcVars)
  
  saveRDS(varsByRepCyc, file=paste(prefix, "Processed/genVars", aof, sep=""))
  
  rm(list=setdiff(ls(), c("prefix", "allOutFiles"))); gc()
}

gvs <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F600NORRANStandard_SelfAllOut.rds")
gvsm <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F600NORRANStandard_MeanSelfAllOut.rds")

#gvs <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F1500NORRANFounderSize1500_SelfAllOut.rds")
#gvsm <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F1500NORRANFounderSize1500_MeanSelfAllOut.rds")

gvs <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F200NORRANFounderSize200_SelfAllOut.rds")
gvsm <- readRDS("Processed/genVars_purgeG2000H1000M1200Q300N500F200NORRANFounderSize200_MeanSelfAllOut.rds")

round(rowMeans(gvs), 5)
round(rowMeans(gvsm), 5)
round(apply(gvs, 1, sd)/sqrt(24), 5)
round(apply(gvsm, 1, sd)/sqrt(24), 5)
sapply(1:6, function(rn) unlist(t.test(x=gvs[rn,], y=gvsm[rn,], paired=T)[1:7]))
