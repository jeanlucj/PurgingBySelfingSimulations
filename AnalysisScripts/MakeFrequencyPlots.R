setwd("~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/PurgingSimulationScript2019")
infiles <- list.files(path="~/Documents/BioHPC/jj332/PurgeDataScripts", pattern="allFreq")
inSuff <- substring(infiles, 8)
fileSuff <- inSuff[1]

# Loop would start here
allFreq <- readRDS(paste0("~/Documents/BioHPC/jj332/PurgeDataScripts/allFreq", fileSuff))

for (fileSuff in inSuff){
  alleleChng <- readRDS(paste0("processed", fileSuff))
  
  # Get only the "main" data
  cycFavR <- t(alleleChng$fAlSt[1, ,])
  cycFavSDR <- t(alleleChng$fAlStSD[1, ,])
  cycDelR <- t(alleleChng$fDel[1, ,])
  cycDelSDR <- t(alleleChng$fDelSD[1, ,])
  
  # Combine data across effect size classes
  effSizes <- colnames(cycFavR)
  effSizes <- sapply(strsplit(effSizes, "_", fixed=T), function(v) v[1])
  effSizes <- as.numeric(substring(effSizes, 3))
  effSzCl <- sort(unique(abs(effSizes)))
  if (length(effSzCl) != 4) stop("There should be 4 size classes")
  effSzName <- c("S", "M", "L", "XL")
  
  procCls <- function(effSz, parmMat, parmMatSD){
    whichCl <- which(abs(effSizes) == effSz)
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
  
  temp <- sapply(effSzCl, procCls, cycFavR, cycFavSDR)
  cycFav <- sapply(1:4*2-1, function(i) temp[[i]])
  cycFavSD <- sapply(1:4*2, function(i) temp[[i]])
  temp <- sapply(effSzCl, procCls, cycDelR, cycDelSDR)
  cycDel <- sapply(1:4*2-1, function(i) temp[[i]])
  cycDelSD <- sapply(1:4*2, function(i) temp[[i]])
  
  cycFavL <- cycFav - cycFavSD/sqrt(24)
  cycFavH <- cycFav + cycFavSD/sqrt(24)
  cycDelL <- cycDel - cycDelSD/sqrt(24)
  cycDelH <- cycDel + cycDelSD/sqrt(24)
  
  stagger <- 0.025
  xlim <- c(0 - stagger, 5 + stagger)
  
  pdf(paste0("pdf", substr(fileSuff, 1, nchar(fileSuff)-4), ".pdf"), width=7, height=6)
  op <- par(mfrow=c(2,2))
  for (locCl in 1:length(effSzCl)){
    ylimDel <- c(min(cycDelL[, locCl]), max(cycDelH[, locCl]))
    ylimFav <- c(min(cycFavL[, locCl]), max(cycFavH[, locCl]))
    par(mar = c(4,4,2,4)) # Leaves some space for the second axis
    plot(0:5 - stagger, cycFav[, locCl], type="l", lwd=2, xlab="Cycle", ylab="Favorable Allele Frequency", col.lab="black", mgp=c(2.3, 1, 0), xlim=xlim, ylim=ylimFav)
    for (cyc in 1:6){
      lines(c(cyc, cyc) - stagger - 1, c(cycFavL[cyc, locCl], cycFavH[cyc, locCl]), lwd=2)
    }
    par(new=T)
    par(xaxt="s")
    plot(0:5 + stagger, cycDel[, locCl], type="l", lwd=2, axes=F, col="dark red", xlab="", ylab="", col.lab="dark red", mgp=c(2.7, 1, 0), xlim=xlim, ylim=ylimDel, main=paste0("Effect Size ", effSzName[locCl]))
    for (cyc in 1:6){
      lines(c(cyc, cyc) + stagger - 1, c(cycDelL[cyc, locCl], cycDelH[cyc, locCl]), lwd=2, col="dark red")
    }
    axis(side=4, col.axis="dark red") # Adds secondary axis
    mtext("Deleterious Homozygote Freq.", side=4, line=2, col="dark red")
  }
  par(op)
  dev.off()
}
