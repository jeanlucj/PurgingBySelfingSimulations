# Make plots 
prePrefix <- "~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/"

########################
# Standard simulations
########################

# Assemble the data
prefix <- "PurgingSimulationScript2019/Schemes/Outputs/"
SelfMean <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANStandard_MeanSelfPlotOut.rds"))$plotData
Self <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANStandard_SelfPlotOut.rds"))$plotData
NoSelf <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANStandard_NoSelfPlotOut.rds"))$plotData
schemeMeans <- cbind(popID=0:5, selfMean=subset(SelfMean, group == 25)$g, self=subset(Self, group == 25)$g, noSelf=subset(NoSelf, group == 25)$g)
SelfMean <- subset(SelfMean, group < 25)
Self <- subset(Self, group < 25)
NoSelf <- subset(NoSelf, group < 25)

# Put selfMean on the bottom, that is, in first position.
# Will require changing orders in legend relative to brkdn.plot col parameter
SelfMean$scheme <- "S1" # make new column
Self$scheme <- "S2"
NoSelf$scheme <- "S3"
Selfmean_Self_NoSelf <- rbind(SelfMean, Self, NoSelf)
Selfmean_Self_NoSelf$group <- paste0("Rep", Selfmean_Self_NoSelf$group)

table(Selfmean_Self_NoSelf$scheme)
head(Selfmean_Self_NoSelf)

# Make the plot Standard
prefix <- "Manuscript/MakeFigures/"
library(plotrix)
mdfunc<-function(x, na.rm=T){return(sd(x, na.rm=T) / sqrt(sum(!is.na(x))))
}

pdf(paste0(prePrefix, prefix, "responseToStandard.pdf"))
par(mar = c(4,4,3,2))
bp <- brkdn.plot(g ~ scheme + popID, data=Selfmean_Self_NoSelf, stagger=0.01, ylim=c(-0.1, 1.2), lwd=2, cex=1, col=c("blue", "orange", "black"), md="mdfunc", xlab="Generation", ylab="Genetic Improvement", main="Response to Standard Schemes", cex.lab=1.3, cex.axis=1.3)
# legend doesn't follow same order as col in
legend("topleft", legend=c("Self", "Self+", "NoSelf"), col=c("orange", "blue","black" ), lwd=2, horiz = T)
dev.off()

########################
# 200 and 1500 Founders
########################
# Assemble the data for F200
prefix <- "PurgingSimulationScript2019/FounderSizes/Outputs/"
SelfMean200 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F200NORRANFounderSize200_MeanSelfPlotOut.rds"))$plotData
Self200 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F200NORRANFounderSize200_SelfPlotOut.rds"))$plotData
NoSelf200 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F200NORRANFounderSize200_NoSelfPlotOut.rds"))$plotData

schemeMeans200 <- cbind(popID=0:5, selfMean=subset(SelfMean200, group == 25)$g, self=subset(Self200, group == 25)$g, noSelf=subset(NoSelf200, group == 25)$g)
SelfMean200 <- subset(SelfMean200, group < 25)
Self200 <- subset(Self200, group < 25)
NoSelf200 <- subset(NoSelf200, group < 25)

SelfMean200$scheme <- "S1"
Self200$scheme <- "S2"
NoSelf200$scheme <- "S3" 

FounderSize200<- rbind(SelfMean200, Self200, NoSelf200)

FounderSize200$Founder <- "F200"
head(FounderSize200)

# Assemble the data for F1500
prefix <- "PurgingSimulationScript2019/FounderSizes/Outputs/"
SelfMean1500 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F1500NORRANFounderSize1500_MeanSelfPlotOut.rds"))$plotData
Self1500 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F1500NORRANFounderSize1500_SelfPlotOut.rds"))$plotData
NoSelf1500 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F1500NORRANFounderSize1500_NoSelfPlotOut.rds"))$plotData

schemeMeans1500 <- cbind(popID=0:5, selfMean=subset(SelfMean1500, group == 25)$g, self=subset(Self1500, group == 25)$g, noSelf=subset(NoSelf1500, group == 25)$g)
SelfMean1500 <- subset(SelfMean1500, group < 25)
Self1500 <- subset(Self1500, group < 25)
NoSelf1500 <- subset(NoSelf1500, group < 25)

SelfMean1500$scheme <- "S1"
Self1500$scheme <- "S2"
NoSelf1500$scheme <- "S3" 

FounderSize1500<- rbind(SelfMean1500, Self1500, NoSelf1500)

FounderSize1500$Founder <- "F1500"
head(FounderSize1500)

# Make the plot F200 F1500
prefix <- "Manuscript/MakeFigures/"
ylim <- range(c(schemeMeans200[,-1], schemeMeans1500[,-1])) # Extreme of the means
ylim <- ylim + c(-0.1, 0.1) # Make space for the error bars

pdf(paste0(prePrefix, prefix, "responseToF200F1500.pdf"), width=14)
par(mar = c(4,4,3,2), mfcol=c(1, 2))
bp <- brkdn.plot(g ~ scheme + popID, data=FounderSize200, stagger=0.01, ylim=ylim, lwd=2, cex=1, col=c("blue", "orange", "black"), md="mdfunc", xlab="Generation", ylab="Genetic Improvement", main="Founder Population Size = 200", cex.lab=1.3, cex.axis=1.3)
# legend doesn't follow same order as col in
legend("topleft", legend=c("Self", "Self+", "NoSelf"), col=c("orange", "blue","black" ), lwd=2, horiz = T)

bp <- brkdn.plot(g ~ scheme + popID, data=FounderSize1500, stagger=0.01, ylim=ylim, lwd=2, cex=1, col=c("blue", "orange", "black"), md="mdfunc", xlab="Generation", ylab="Genetic Improvement", main="Founder Population Size = 1500", cex.lab=1.3, cex.axis=1.3)
# legend doesn't follow same order as col in
legend("topleft", legend=c("Self", "Self+", "NoSelf"), col=c("orange", "blue","black" ), lwd=2, horiz = T)
dev.off()

########################
# Additive gene action
########################
# Assemble the data: Additive
prefix <- "PurgingSimulationScript2019/AdditiveModel/Outputs/"
AdditiveF200 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F200NORRANAdditiveModelFounder200SelfPlotOut.rds"))$plotData
AdditiveF600 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANAdditiveModelFounder600SelfPlotOut.rds"))$plotData
AdditiveF1500 <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F1500NORRANAdditiveModelFounder1500SelfPlotOut.rds"))$plotData

AllAdditiveMM3 <- subset(AdditiveF200, group!=25)
AllAdditiveMM2 <- subset(AdditiveF600, group!=25)
AllAdditiveMM1 <- subset(AdditiveF1500, group!=25)

AllAdditiveMM3$Founder <- rep("C0200",nrow(AllAdditiveMM3))
AllAdditiveMM2$Founder <- rep("A0600",nrow(AllAdditiveMM2))
AllAdditiveMM1$Founder <- rep("B1500",nrow(AllAdditiveMM1))

AllADDFounder <- rbind(AllAdditiveMM1, AllAdditiveMM2,AllAdditiveMM3)

# Make the plot F200 F1500
prefix <- "Manuscript/MakeFigures/"
# plotted in order 1500, 600, 200.
pdf(paste0(prePrefix, prefix, "responseToAdditive.pdf"))
par(mar = c(4,4,3,2))
bp <- brkdn.plot(g ~ Founder + popID, data=AllADDFounder, stagger=0.01, lwd=2, cex=1, col=c("blue", "brown", "orange"), md="mdfunc", xlab="Generation", ylab="Genetic Improvement", main="Additive Model by Founder Population Size", cex.lab=1.3, cex.axis=1.3)
legend("topleft",legend=c("200","600", "1500"), col=c("orange", "blue", "brown"), lwd = 2, horiz = T)
dev.off()

########################
# High selection intensity
########################
# Assemble the data: Additive
prefix <- "PurgingSimulationScript2019/SelectionIntensity/Outputs/"
HselectionIntenSM <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANHighSelectionIntensityF600_MeanSelfPlotOut.rds"))$plotData
HselectionIntenS <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANHighSelectionIntensityF600_SelfPlotOut.rds"))$plotData
HselectionIntenNS <- readRDS(paste0(prePrefix, prefix, "_purgeG2000H1000M1200Q300N500F600NORRANHighSelectionIntensityF600_NoSelfPlotOut.rds"))$plotData

schemeMeansHigh <- cbind(popID=0:5, selfMean=subset(HselectionIntenSM, group == 25)$g, self=subset(HselectionIntenS, group == 25)$g, noSelf=subset(HselectionIntenNS, group == 25)$g)

HselectionIntenSM$scheme <- "S1"
HselectionIntenS$scheme <- "S2"
HselectionIntenNS$scheme <- "S3" 

AllHighSelectInten<- rbind(HselectionIntenSM, HselectionIntenS, HselectionIntenNS)
AllHighSelectInten$Selection<- "HighIntensity"
AllHighSelectIntenNoMean <- subset(AllHighSelectInten, group!=25)

# Make the plot Standard
prefix <- "Manuscript/MakeFigures/"

pdf(paste0(prePrefix, prefix, "responseAtHighSelectionIntensity.pdf"))
par(mar = c(4,4,3,2))
bp <- brkdn.plot(g ~ scheme + popID, data=AllHighSelectIntenNoMean, stagger=0.01, lwd=2, cex=1, col=c("blue", "orange", "black"), md="mdfunc", xlab="Generation", ylab="Genetic Improvement", main="Response at High Selection Intensity", cex.lab=1.3, cex.axis=1.3)
# legend doesn't follow same order as col in
legend("topright", legend=c("Self", "Self+", "NoSelf"), col=c("orange", "blue","black" ), lwd=2, horiz = T)
dev.off()

########################
# Pattern of no progress after the first selection
# Homozygote deleterious frequency increases even as does the favorable allele
########################
# WARNING: This was messed up by renaming the effect sizes...
# Assemble the data
inFiles <- c("processed_F200NORRANFounderSize200_SelfAllOut.rds", "processed_F600NORRANStandard_SelfAllOut.rds", "processed_F1500NORRANFounderSize1500_SelfAllOut.rds", "processed_F600NORRANAdditiveModelFounder600SelfAllOut.rds", "processed_F600NORRANHighSelectionIntensityF600_SelfAllOut.rds")

preFix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/Processed/"
makeFreqChngPlot <- function(fileName){
  print(fileName)
  alleleChng <- readRDS(paste0(preFix, fileName))
  
  # Get only the "main" data
  cycFavR <- t(alleleChng$fAlSt[1, ,]) # Frequency of the favorable allele
  cycFavSDR <- t(alleleChng$fAlStSD[1, ,])
  cycDelR <- t(alleleChng$fDel[1, ,]) # Frequency of the deleterious homozygote
  cycDelSDR <- t(alleleChng$fDelSD[1, ,])
  
  # Combine data across effect size classes
  effSizes <- colnames(cycFavR)
  # Bad work around.  Should redo the _processed file.
  if (length(grep("L", effSizes)) == 0){
    effSizes <- c("efM_at0.1", "efL_at0.05", "efX_at0.02", "efS_at0.9", "efM_at0.9", "efM_at0.95")
    colnames(cycFavR) <- effSizes
    colnames(cycFavSDR) <- effSizes
    colnames(cycDelR) <- effSizes
    colnames(cycDelSDR) <- effSizes
  } 
  effSizes <- sapply(strsplit(effSizes, "_", fixed=T), function(v) v[1])
  effSizes <- substring(effSizes, 3, 3)
  # effSizes[effSizes == "X"] <- "XL"
  # effSzOrder <- 1:4
  namesEffBySz <- c("S", "M", "L", "X")
  # effSzCl <- sort(unique(abs(effSizes)))
  # if (length(effSzCl) != 4) stop("There should be 4 size classes")
  # effSzName <- c("S", "M", "L", "XL")
  
  procCls <- function(effNm, parmMat, parmMatSD){
    whichCl <- grep(effNm, colnames(parmMat))
    if (length(whichCl) > 1){
      combineEff <- function(valMat, sdMat){
        sdMat <- 1 / (sdMat^2)
        cOneRow <- function(i){
          return((valMat[i,] %*% sdMat[i,]) / sum(sdMat[i,]))
        }
        return(list(parmVec=sapply(1:nrow(valMat), cOneRow), parmVecSD=apply(sdMat, 1, function(v) sqrt(1 / sum(v)))))
      }
      return(combineEff(parmMat[,whichCl], parmMatSD[,whichCl]))
    } else{
      return(list(parmVec=parmMat[,whichCl], parmVecSD=parmMatSD[,whichCl]))
    }
  }#END procCls
  
  temp <- sapply(namesEffBySz, procCls, cycFavR, cycFavSDR)
  cycFav <- sapply(1:4*2-1, function(i) temp[[i]])
  cycFavSD <- sapply(1:4*2, function(i) temp[[i]])
  temp <- sapply(namesEffBySz, procCls, cycDelR, cycDelSDR)
  cycDel <- sapply(1:4*2-1, function(i) temp[[i]])
  cycDelSD <- sapply(1:4*2, function(i) temp[[i]])
  
  cycFavL <- cycFav - cycFavSD/sqrt(24)
  cycFavH <- cycFav + cycFavSD/sqrt(24)
  cycDelL <- cycDel - cycDelSD/sqrt(24)
  cycDelH <- cycDel + cycDelSD/sqrt(24)
  
  stagger <- 0.025
  xlim <- c(0 - stagger, 5 + stagger)
  namesEffBySz <- c("S", "M", "L", "XL")
  
  pdf(paste0(prePrefix, "Manuscript/MakeFigures/freqPlot", substr(fileName, 1, nchar(fileName)-4), ".pdf"), width=7, height=6)
  op <- par(mfrow=c(2,2))
  for (locCl in 1:length(namesEffBySz)){
    ylimDel <- c(min(cycDelL[, locCl]), max(cycDelH[, locCl]))
    ylimFav <- c(min(cycFavL[, locCl]), max(cycFavH[, locCl]))
    par(mar = c(4,4,2,4)) # Leaves some space for the second axis
    plot(0:5 - stagger, cycFav[, locCl], type="l", lwd=2, xlab="Cycle", ylab="Favorable Allele Frequency", col.lab="black", mgp=c(2.3, 1, 0), xlim=xlim, ylim=ylimFav)
    for (cyc in 1:6){
      lines(c(cyc, cyc) - stagger - 1, c(cycFavL[cyc, locCl], cycFavH[cyc, locCl]), lwd=2)
    }
    par(new=T)
    par(xaxt="s")
    plot(0:5 + stagger, cycDel[, locCl], type="l", lwd=2, axes=F, col="dark red", xlab="", ylab="", col.lab="dark red", mgp=c(2.7, 1, 0), xlim=xlim, ylim=ylimDel, main=paste0("Effect Size ", namesEffBySz[locCl]))
    for (cyc in 1:6){
      lines(c(cyc, cyc) + stagger - 1, c(cycDelL[cyc, locCl], cycDelH[cyc, locCl]), lwd=2, col="dark red")
    }
    axis(side=4, col.axis="dark red") # Adds secondary axis
    mtext("Deleterious Homozygote Freq.", side=4, line=2, col="dark red")
  }
  par(op)
  dev.off()
}
dummy <- sapply(inFiles, makeFreqChngPlot)

########################
# Pattern that purge selection only helps in the first selection
# The two hypotheses are
# 1. The ability of selfing to increase favorable allele frequencies fails
#    after the first selection
# 2. The main selection drops in efficacy after the first selection under
#    selfing but not under no selfing
# The latter ends up being true. I need the changes in allele frequencies
# under each scheme after main and purge selection events
########################
setwd("~/Documents/BioHPC/jj332/PurgeDataScripts/Processed")
resPV6 <- readRDS(file="AverageShiftWgtFavAlFreqOverFounderPops.rds")
names(resPV6) <- c("NoSelf", "Self")
colnames(resPV6$NoSelf) <- colnames(resPV6$Self) <- paste0("Sel", 1:5)
rownames(resPV6$NoSelf) <- rownames(resPV6$Self) <- c("Main", "Purge", "MainSD", "PurgeSD")
mainNoSelf <- resPV6$NoSelf["Main",]
mainNoSelfSE <- resPV6$NoSelf["MainSD",]/sqrt(24)
mainSelf <- resPV6$Self["Main",]
mainSelfSE <- resPV6$Self["MainSD",]/sqrt(24)
purgeSelf <- resPV6$Self["Purge",]
purgeSelfSE <- resPV6$Self["PurgeSD",]/sqrt(24)
mainSE <- (mainNoSelfSE + mainSelfSE) / 2

# Lefthand scale is going to be main change in allele frequency for Self and NoSelf
# Righthand scale is going to be purge change in allele frequency for Self
pdf(paste0(prePrefix, "Manuscript/MakeFigures/afterSel1.pdf"))
stagger <- 0.025
xlim <- c(1 - stagger, 5 + stagger)
# ylim <- c(min(c(mainNoSelf - mainSE, mainSelf - mainSE, purgeSelf - purgeSelfSE)), max(c(mainNoSelf + mainSE, mainSelf + mainSE, purgeSelf + purgeSelfSE)))
ylim <- c(0, max(c(mainNoSelf + mainSE, mainSelf + mainSE)))
par(mar = c(4,4,2,4)) # Leaves some space for the second axis
# Lines for NoSelf -- Main
plot(1:5 + stagger, mainNoSelf, type="l", lwd=2, xlab="Selection Event", ylab="Main Allele Freq. Change", col.lab="black", mgp=c(2.3, 1, 0), xlim=xlim, ylim=ylim, lty=2, cex.lab=1.3, cex.axis=1.3)
for (selEv in 1:5){
  lines(c(selEv, selEv) + stagger, mainNoSelf[selEv]+c(-mainSE[selEv], mainSE[selEv]), lwd=2)
}
# Lines for Self -- Main
lines(1:5 - stagger, mainSelf, lwd=2)
for (selEv in 1:5){
  lines(c(selEv, selEv) - stagger, mainSelf[selEv]+c(-mainSE[selEv], mainSE[selEv]), lwd=2)
}

# Lines for Purge
par(new=T)
par(xaxt="s")
ylimPur <- ylim / 2
plot(1:5 - stagger, purgeSelf, type="l", lwd=2, axes=F, col="dark red", xlab="", ylab="", col.lab="dark red", mgp=c(2.7, 1, 0), xlim=xlim, ylim=ylimPur)
for (selEv in 1:5){
  lines(c(selEv, selEv) - stagger, purgeSelf[selEv]+c(-purgeSelfSE[selEv], purgeSelfSE[selEv]), lwd=2, col="dark red")
}
axis(side=4, col.axis="dark red", cex.axis=1.3) # Adds secondary axis
mtext("Purge Allele Freq. Change", side=4, line=2.2, col="dark red", cex=1.3)
dev.off()

########################
# Purge selection does not increase the favorable allele freq.
# as much as we might expect based on the increase in genotypic value
########################
setwd("~/Documents/BioHPC/jj332/PurgeDataScripts/Processed")
procFiles <- list.files(pattern="processed_")
procFiles <- c("processed_F1500NORRANAdditiveModelFounder1500SelfAllOut.rds", "processed_F1500NORRANFounderSize1500_MeanSelfAllOut.rds", "processed_F1500NORRANFounderSize1500_SelfAllOut.rds", "processed_F200NORRANAdditiveModelFounder200SelfAllOut.rds", "processed_F200NORRANFounderSize200_MeanSelfAllOut.rds", "processed_F200NORRANFounderSize200_SelfAllOut.rds", "processed_F600NORRANAdditiveModelFounder600SelfAllOut.rds", "processed_F600NORRANStandard_MeanSelfAllOut.rds", "processed_F600NORRANStandard_SelfAllOut.rds")
parmsToProc <- c("fChng", "fLost")
parmsToProc <- c("fChng") # Use this for resPV6

# !!!WARNING!!! The wgts should be in alphabetical order of the effect sizes
# These effect sizes are now S, M, L, X.  So that order goes L, M, M, M, S, X
# ALSO, the additive now only has FOUR effects (collapse the three M=2)
# !!!WARNING!!! This code used for both "MainPurgeAlleleFreqChange.pdf" and for
# figure showing decreased efficacy of main selection under self
# The two depend on what gets returned
procProc <- function(fileName){
  isAdd <- length(grep("Additive", fileName)) > 0
  if (isAdd) wgt <- c(4, 2, 1, 8) else wgt <- c(4, 2, 2, 2, 1, 8)
  allParms <- readRDS(fileName)
  doParm <- function(whichParm){
    parmMat <- allParms[[whichParm]]
    parmMatSD <- allParms[[paste0(whichParm, "SD")]]
    parmMat <- parmMat[,,-6] # Take off the last cycle: no change or loss
    parmMatSD <- parmMatSD[,,-6]
    parmMatByCyc <- apply(parmMat, c(1, 3), weighted.mean, w=wgt)
    varMat <- parmMatSD^2
    varMatByCyc <- apply(varMat, c(1, 3), weighted.mean, w=wgt^2)*sum(wgt^2)/sum(wgt)^2
    meanAcrossCycles <- function(rn){ # Sometimes the variance is zero
      if (any(varMatByCyc[rn,]==0)) return(mean(parmMatByCyc[rn,]))
      else return(weighted.mean(parmMatByCyc[rn,], w=1/varMatByCyc[rn,]))
    }
    mp <- sapply(1:2, meanAcrossCycles)
    mainSD <- sqrt(1 / sum(1/varMatByCyc[1,]))
    purgeSD <- sqrt(1 / sum(1/varMatByCyc[2,]))
    # Choose return depending on whether you want averaged over cycles or not
    return(c(main=mp[1], purge=mp[2], mainSD=mainSD, purgeSD=purgeSD))
    # return(rbind(parmMatByCyc, varMatByCyc)) # Use this for resPV6
  }
  return(unlist(lapply(parmsToProc, doParm)))
  # return(lapply(parmsToProc, doParm)) # Use this for resPV6
}
pV <- sapply(procFiles, procProc)
colNoSelf <- grep("_NoSelf", colnames(pV))
procAlFreq <- pV
if (length(colNoSelf) > 0) procAlFreq <- pV[,-colNoSelf]

setwd("~/Documents/BioHPC/jj332/PurgeDataScripts/MainPurge")
dataFiles <- list.files(pattern="data_")
procValsFile <- function(fileName){
  allVals <- readRDS(fileName)
  getMeanSD <- function(valsDF){
    valsDF$group <- as.numeric(valsDF$group)
    valsDF$popID <- as.numeric(valsDF$popID)
    valsDF <- valsDF[valsDF$group < 25 & valsDF$popID > 0,]
    gDiff <- diff(valsDF$g)[1:120 * 2 - 1] # lagged difference
    replication <- as.factor(rep(1:24, each=5)) # 24 reps of 5 cycles
    byRep <- tapply(gDiff, replication, mean)
    return(c(mean=mean(byRep), sd=sd(byRep)))
  }
  main <- getMeanSD(allVals$mainVals)
  purge <- getMeanSD(allVals$purgeVals)
  return(c(main=main[1], purge=purge[1], mainSD=main[2], purgeSD=purge[2]))
}
pV <- sapply(dataFiles, procValsFile)
procGenVal <- pV[,-grep("_NoSelf", colnames(pV))]

# Going to try looking at this differently, with genotypic value changes
# on one axis and allele frequency changes on the other
procAlFr <- procAlFreq[1:4,]
# reshape procAlFr and procGenVal
hasVal <- function(val, vec){
  hs <- grep(val, vec)
  res <- rep(F, length(vec))
  res[hs] <- T
  return(res)
}
reshapeProc <- function(procDF){
  tp <- t(procDF)
  hasAdd <- hasVal("Additive", colnames(procDF))
  hasSlf <- hasVal("_Self", colnames(procDF))
  hasMSl <- hasVal("_MeanSelf", colnames(procDF))
  genMod <- c("Additive", "Self", "MeanSelf")[hasAdd+hasSlf*2+hasMSl*3]
  f200 <- hasVal("F200", colnames(procDF))
  f600 <- hasVal("F600", colnames(procDF))
  f1500 <- hasVal("F1500", colnames(procDF))
  foundSz <- c("F200", "F600", "F1500")[f200+f600*2+f1500*3]
  mainPurge <- rep(c("main", "purge"), each=nrow(tp))
  tp <- rbind(tp[,c(1,3)], tp[,c(2,4)])
  res <- data.frame(mainPurge, genMod, foundSz, tp)
  colnames(res)[4:5] <- c("chng", "chngSD")
  return(res)
}

resAlFr <- reshapeProc(procAlFr)
resGenVal <- reshapeProc(procGenVal)

plotCol <- as.numeric(resAlFr$genMod)
plotPch <- 14 + as.numeric(resAlFr$foundSz)
xSE <- resGenVal$chngSD / sqrt(24); ySE <- resAlFr$chngSD / sqrt(24)
xRange <- rbind(resGenVal$chng - xSE, resGenVal$chng + xSE)
yRange <- rbind(resAlFr$chng - ySE, resAlFr$chng + ySE)
xlim <- range(c(xRange)); xlim[1] <- 0
ylim <- range(c(yRange)); ylim[1] <- 0
xlab <- "Genotypic Value Change"
ylab <- "Favorable Allele Freq. Change"
pdf(paste0(prePrefix, "Manuscript/MakeFigures/AlleleFreqVsGenoVal.pdf"))
plot(resGenVal$chng, resAlFr$chng, col=plotCol, pch=plotPch, xlab=xlab, ylab=ylab, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlim=xlim, ylim=ylim)
# Do the stdErr
for (pt in 1:18){
  lines(xRange[,pt], rep(resAlFr$chng[pt], 2), lwd=2, col=plotCol[pt])
  lines(rep(resGenVal$chng[pt], 2), yRange[,pt], lwd=2, col=plotCol[pt])
}
dev.off()

########################
# Population "rescue" using the Self+ scheme
# Homozygote deleterious frequency increases even as does the favorable allele
########################
# WARNING: This was messed up by renaming the effect sizes...
# Assemble the data
inFiles <- list(c("processed_F600NORRANHighSelectionIntensityF600_SelfAllOut.rds", "processed_F600NORRANHighSelectionIntensityF600_MeanSelfAllOut.rds"), c("processed_F200NORRANFounderSize200_SelfAllOut.rds", "processed_F200NORRANFounderSize200_MeanSelfAllOut.rds"), c("processed_F1500NORRANFounderSize1500_MeanSelfAllOut.rds", "processed_F1500NORRANFounderSize1500_SelfAllOut.rds"))

titles <- c("High Intensity Self", "High Intensity Self+", "Founder Size 200 Self", "Founder Size 200 Self+", "Founder Size 1500 Self", "Founder Size 1500 Self+")
titleIdx <- 1

preFix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/Processed/"
makeFreqChngPlotM <- function(fileNames){
  print(fileNames)
  # limits for _two_ plots
  ylimDel <- c(10, -10)
  ylimFav <- c(10, -10)
  
  makeCoordFromFile <- function(fileName){
    alleleChng <- readRDS(paste0(preFix, fileName))
    
    # Get only the "main" data
    cycFavR <- t(alleleChng$fAlSt[1, ,]) # Frequency of the favorable allele
    cycFavSDR <- t(alleleChng$fAlStSD[1, ,])
    cycDelR <- t(alleleChng$fDel[1, ,]) # Frequency of the deleterious homozygote
    cycDelSDR <- t(alleleChng$fDelSD[1, ,])
    
    effSizes <- colnames(cycFavR)
    # Bad work around for renaming the effect sizes.  Should redo the _processed file.
    if (length(grep("L", effSizes)) == 0){
      effSizes <- c("efM_at0.1", "efL_at0.05", "efX_at0.02", "efS_at0.9", "efM_at0.9", "efM_at0.95")
      colnames(cycFavR) <- effSizes
      colnames(cycFavSDR) <- effSizes
      colnames(cycDelR) <- effSizes
      colnames(cycDelSDR) <- effSizes
    } 
    effSizes <- sapply(strsplit(effSizes, "_", fixed=T), function(v) v[1])
    effSizes <- substring(effSizes, 3, 3)
    namesEffBySz <- c("S", "M", "L", "X")
    
    # function to combine the different versions of the M effect size class
    procCls <- function(effNm, parmMat, parmMatSD){
      whichCl <- grep(effNm, colnames(parmMat))
      if (length(whichCl) > 1){
        combineEff <- function(valMat, sdMat){
          sdMat <- 1 / (sdMat^2)
          cOneRow <- function(i){
            return((valMat[i,] %*% sdMat[i,]) / sum(sdMat[i,]))
          }
          return(list(parmVec=sapply(1:nrow(valMat), cOneRow), parmVecSD=apply(sdMat, 1, function(v) sqrt(1 / sum(v)))))
        }
        return(combineEff(parmMat[,whichCl], parmMatSD[,whichCl]))
      } else{
        return(list(parmVec=parmMat[,whichCl], parmVecSD=parmMatSD[,whichCl]))
      }
    }#END procCls
    
    temp <- sapply(namesEffBySz, procCls, cycFavR, cycFavSDR)
    cycFav <- sapply(1:4*2-1, function(i) temp[[i]])
    cycFavSD <- sapply(1:4*2, function(i) temp[[i]])
    temp <- sapply(namesEffBySz, procCls, cycDelR, cycDelSDR)
    cycDel <- sapply(1:4*2-1, function(i) temp[[i]])
    cycDelSD <- sapply(1:4*2, function(i) temp[[i]])
    
    cycFavL <- cycFav - cycFavSD/sqrt(24)
    cycFavH <- cycFav + cycFavSD/sqrt(24)
    cycDelL <- cycDel - cycDelSD/sqrt(24)
    cycDelH <- cycDel + cycDelSD/sqrt(24)
    
    locCl <- 2 # For the M effect size
    ylimFav <<- c(min(c(ylimFav[1], cycFavL[, locCl])), max(c(ylimFav[2], cycFavH[, locCl])))
    ylimDel <<- c(min(c(ylimDel[1], cycDelL[, locCl])), max(c(ylimDel[2], cycDelH[, locCl])))
    
    return(list(cycFav, cycDel, cycFavL, cycFavH, cycDelL, cycDelH))
  }#END makeCoord
  
  coordLists <- lapply(fileNames, makeCoordFromFile)
  
  makePlotFromCoord <- function(coordList){
    stagger <- 0.025
    xlim <- c(0 - stagger, 5 + stagger)
    namesEffBySz <- c("S", "M", "L", "XL")
    locCl <- 2 # For the M effect size
    
    cycFav <- coordList[[1]]; cycDel <- coordList[[2]]
    cycFavL <- coordList[[3]]; cycFavH <- coordList[[4]]
    cycDelL <- coordList[[5]]; cycDelH <- coordList[[6]]
    
    plot(0:5 - stagger, cycFav[, locCl], type="l", lwd=2, xlab="Cycle", ylab="Favorable Allele Frequency", col.lab="black", mgp=c(2.3, 1, 0), xlim=xlim, ylim=ylimFav, cex.axis=1.3, cex.lab=1.3)
    for (cyc in 1:6){
      lines(c(cyc, cyc) - stagger - 1, c(cycFavL[cyc, locCl], cycFavH[cyc, locCl]), lwd=2)
    }
    par(new=T)
    par(xaxt="s")
    plot(0:5 + stagger, cycDel[, locCl], type="l", lwd=2, axes=F, col="dark red", xlab="", ylab="", col.lab="dark red", mgp=c(2.7, 1, 0), xlim=xlim, ylim=ylimDel, main=titles[titleIdx])
    titleIdx <<- titleIdx + 1
    for (cyc in 1:6){
      lines(c(cyc, cyc) + stagger - 1, c(cycDelL[cyc, locCl], cycDelH[cyc, locCl]), lwd=2, col="dark red")
    }
    axis(side=4, col.axis="dark red", cex.axis=1.3) # Adds secondary axis
    mtext("Deleterious Homozygote Freq.", side=4, line=2, col="dark red", cex=0.8)
    return(NULL)
  }#END makePlot
  
  dummy <- lapply(coordLists, makePlotFromCoord)
}

pdf(paste0(prePrefix, "Manuscript/MakeFigures/SelfPlusRescue.pdf"), width=7, height=9)
op <- par(mfrow=c(3,2))
par(mar = c(4,4,2,4)) # Leaves some space for the second axis

dummy <- sapply(inFiles, makeFreqChngPlotM)

par(op)
dev.off()

