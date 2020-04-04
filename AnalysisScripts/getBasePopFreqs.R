###############################################################################
# Scripts to look at the underlying mechanisms in terms of allele and genotype
# frequencies for loci across individuals of for individuals across classes
# of loci
###############################################################################

# Set correct folders depending on compute environment
bioHPC <- TRUE
if (bioHPC){
  setwd("/home/jj332/PurgeDataScripts/")
  prefix <- "/home/jj332/PurgeDataScripts/"
} else{
  setwd("~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/PurgingSimulationScript2019")
  prefix <- "~/Documents/BioHPC/jj332/PurgeDataScripts/"
}

library(lme4)
# stdSelf <- readRDS(paste0(prefix, "AllOut/_purgeG2000H1000M1200Q300N500F600NORRANStandard_SelfAllOut.rds"))
# Number of favorable versus deleterious fixed
# Degree of dominance of what is still segregating
# Function to handle one simulation replication at a time
# Possible values for analysisToDo are
# locusClass. A locusClass is the combination of effect and degree of dominance
#     Analysis gives the means for each base population across all individuals
#     and loci for the favorable allele frequency for each locusClass
# freqDat. Here you want frequency data by locus. You create a matrix with loci
#     in rows and start and ending frequencies of alleles and genotypes, per
#     cycle, and whether it's a main or a purge cycle
# selCoef. Some kind of regression on pq that approximates the selection
#     coefficient, also locus by locus
# byInd. Frequency data by locusClass, within individuals. So the matrix is
#     with individuals in rows, four columns per locusClass, with the genotype
#     and allele frequencies. From there, calculate the weighted frequencies
calcLocParms <- function(bslData, analysisToDo="locusClass"){
  genoBSL <- bslData$geno[,bslData$mapData$effectivePos]
  # The output I want is a table showing the frequencies of the alleles
  # classified by effect and actionType
  # A locus class is the combo of effect and actionType
  makeLocClass <- function(effects, actionType){
    effCl <- sort(unique(effects))
    actCl <- sort(unique(actionType))
    locCl <- list()
    for (eff in effCl){
      for (act in actCl){
        inCl <- which(effects == eff & actionType == act)
        if (length(inCl) > 0) locCl <- c(locCl, list(list(whichInClass=inCl, effect=eff, actionType=act)))
      }
    }
    return(locCl)
  }
  # There are six classes
  locClasses <- makeLocClass(bslData$mapData$effects, bslData$mapData$actionType)
  
  # The name includes the effect and actionType
  makeLocClName <- function(locCl){
    effCl <- c("S", "M", "L", "X")[(abs(locCl$effect) > 0.18) + (abs(locCl$effect) > 0.09) + (abs(locCl$effect) > 0.045) + 1]
    paste0("ef", effCl, "_at", round(locCl$actionType, 2))
  }
  locClNames <<- sapply(locClasses, makeLocClName)
  nClasses <<- length(locClNames)
  
  if (analysisToDo == "byInd"){
    # Calculate the frequency of genotypes at in classes for each individual
    # This corrects so that the favorable allele is always +1
    calcGenoFreqByInd <- function(locClass, geno){
      loc <- locClass$whichInClass
      nInd <- nrow(geno) / 2
      dosages <- sign(locClass$effect) * 
        (geno[1:nInd * 2, loc] + geno[1:nInd * 2 - 1, loc]) / 2 + 1
      genoFreq <- function(vec){
        counts <- table(vec)
        counts <- counts / sum(counts)
        freq <- counts %*% as.integer(names(counts))/2
        missGeno <- setdiff(c("0", "1", "2"), names(counts))
        toRet <- integer(length(missGeno))
        names(toRet) <- missGeno
        toRet <- c(counts, toRet)
        return(c(toRet[order(names(toRet))], freq=freq))
      }
      return(t(apply(dosages, 1, genoFreq)))
    }
    freqInInd <- sapply(locClasses, calcGenoFreqByInd, geno=genoBSL)
    freqInInd <- matrix(freqInInd, nrow=nrow(genoBSL)/2)
    colnames(freqInInd) <- paste0(rep(locClNames, each=4), c("_0", "_1", "_2", "Al"))
    return(freqInInd)
  }
  
  # Calculate the frequency of genotypes at loci individually in a class
  # This corrects so that the favorable allele is always +1
  # locClass is the locus class created by makeLocClass
  # geno is the matrix of allele states for the whole simulation
  # ind is a vector of individuals representing a population
  # returns a 4 x nLoc matrix where nLoc is the number of loci in the class
  calcGenoFreq <- function(locClass, geno, ind){
    genoFreqByLoc <- function(loc){
      counts <- table(sign(locClass$effect) * (geno[ind * 2, loc] + geno[ind * 2 - 1, loc]))
      names(counts) <- (as.integer(names(counts)) + 2) / 2
      counts <- counts / sum(counts)
      freq <- counts %*% as.integer(names(counts))/2
      missGeno <- setdiff(c("0", "1", "2"), names(counts))
      toRet <- integer(length(missGeno))
      names(toRet) <- missGeno
      toRet <- c(counts, toRet)
      return(c(toRet[order(names(toRet))], freq=freq))
    }
    return(t(sapply(locClass$whichInClass, genoFreqByLoc)))
  }
  
  # Calculate the frequency of loci individually in a class
  # This corrects so that the favorable allele is always +1
  calcLocFreq <- function(locClass, geno, ind){
    gam <- c(ind * 2, ind * 2 - 1)
    loc <- locClass$whichInClass
    # colMeans here: each column is a locus
    return((sign(locClass$effect) * colMeans(geno[gam, loc]) + 1)/2)
  }
  
  # Calculate the frequency of all loci combined in a class
  # This corrects so that the favorable allele is always +1
  calcClassFreq <- function(locClass, geno, ind){
    gam <- c(ind * 2, ind * 2 - 1)
    loc <- locClass$whichInClass
    # mean here: all columns together for a class
    return((sign(locClass$effect) * mean(geno[gam, loc]) + 1)/2)
  }
  
  # get frequencies for each locus or each class of loci
  # Possible values for baseOrSelected are
  #    base: all the individuals made as progeny in the given base cycle
  #    selected: individuals that were selected into a population with popID
  # Possible values for locOrClass
  #    loc: calculate the frequency of alleles locus by locus
  #    geno: calculate the frequency of genotypes locus by locus
  #    class: calculate the frequency of alleles for a locus class
  # A base population is when it is first made as progeny
  # popID can change if the base is selected
  getBasePopFreqs <- function(popID, baseOrSelected="base", locOrClass="loc"){
    if (baseOrSelected == "base"){
      targetIdx <- bslData$genoRec$basePopID
    } else{
      targetIdx <- bslData$genoRec$popID
    }
    ind <- which(targetIdx == popID)
    return(switch(locOrClass,
                  class = sapply(locClasses, calcClassFreq, geno=genoBSL, ind=ind),
                  loc = sapply(locClasses, calcLocFreq, geno=genoBSL, ind=ind),
                  geno = sapply(locClasses, calcGenoFreq, geno=genoBSL, ind=ind))
    )
  }
  
  if (analysisToDo == "locusClass"){
    locClFreq <- t(sapply(sort(unique(bslData$genoRec$basePopID)), getBasePopFreqs, locOrClass="class"))
    colnames(locClFreq) <- locClNames
    return(locClFreq)
  }
  
  # This function takes data from the lists of dataframes baseFreq, seleFreq
  # and puts it all together in one dataframe of a locus class
  makeFreqData <- function(locCl, mainPurge=NULL){
    # Make the dataset
    freqDat <- data.frame()
    
    addToFreqDat <- function(mainPurge){
      className <- locClNames[locCl]
      locCl <- locCl + ifelse(mainPurge == "main", 0, nClasses)
      idx <- locCl + rep(0:4 * (2*nClasses), each=length(locCl)) # Combine selection cycles
      for (i in idx){
        freqDat <<- rbind(freqDat, cbind(data.frame(mainPurge=mainPurge, locusClass=className, cycle=paste0("C", (i-1) %/% (2*nClasses)), stringsAsFactors=F), baseFreq[[i]], seleFreq[[i]]))
      }
      # The last population was not selected, so there is no selected frequency
      if (mainPurge == "main"){
        noSF <- array(NA, dim=dim(baseFreq[[i]]))
        colnames(noSF) <- colnames(seleFreq[[i]])
        i <- i + 2*nClasses
        freqDat <<- rbind(freqDat, cbind(data.frame(mainPurge=mainPurge, locusClass=className, cycle=paste0("C", (i-1) %/% (2*nClasses)), stringsAsFactors=F), baseFreq[[i]], noSF))
      }
    }#END addToFreqDat
    if (is.null(mainPurge)){
      for (mainPurge in c("main", "purge")) addToFreqDat(mainPurge)
    } else addToFreqDat(mainPurge)
    return(freqDat)
  }#END makeFeqDat
  
  ###########################################
  # This is the work horse
  ###########################################
  if (analysisToDo == "freqDat"){
    # Lists: six loc classes x number of populations
    # Loc classes cycle then the populations
    # Loc classes in order: 
    # efX_at0.02 efL_at0.05 efM_at0.1
    # efS_at0.9 efM_at0.9 efM_at0.95
    # Once class 3 is flipped, it resembles 5 and 6
    # Frequencies at each locus in each locus class across base populations
    baseFreq <- sapply(sort(unique(bslData$genoRec$basePopID)), getBasePopFreqs, locOrClass="geno")
    # Frequencies at each locus in each locus class across selected populations
    seleFreq <- sapply(c(1, 3, 5, 7, 9, 11, 14, 16, 20, 22), getBasePopFreqs, baseOrSelected="selected", locOrClass="geno")
    
    # Combine allele frequency data for a locus class across cycles
    # locCl is an integer vector in the range of 1 to 6 indexing the loc class
    
    allFreqDat <- lapply(1:nClasses, makeFreqData)
    return(allFreqDat)
  }#END analysisToDo == "freqDat"
  
  # Did not use this in the end
  if (analysisToDo == "selCoef"){
    # Frequencies at each locus in each locus class across base populations
    baseFreq <- sapply(sort(unique(bslData$genoRec$basePopID)), getBasePopFreqs, locOrClass="loc")
    # Frequencies at each locus in each locus class across selected populations
    seleFreq <- sapply(c(1, 3, 5, 7, 9, 11, 14, 16, 20, 22), getBasePopFreqs, baseOrSelected="selected", locOrClass="loc")
    
    # Calculate the selection coefficient by regressing on pq
    # That forces the regression curve to zero at p=0 and p=1
    # It follows more or less the expected change in allele frequency F&M p. 28
    calcSelCoef <- function(locCl, mainPurge, doPlot=TRUE){
      freqDat <- makeFreqData(locCl, mainPurge)
      colnames(freqDat)[4:5] <- c("end", "start")
      deltaP <- freqDat$end - freqDat$start
      startPQ <- freqDat$start*(1 - freqDat$start)
      fitDP <- lm(deltaP ~ -1 + startPQ)
      sVal <- coef(fitDP)
      seSval <- summary(fitDP)$coefficients[2]
      if (doPlot){
        plot(freqDat$start, deltaP, cex=0.8, pch=16, col=as.integer(substr(freqDat$cycle, start=2, stop=2)), main=paste(mainPurge, locCl))
        p <- 0:100 / 100
        dp <- sVal * p * (1 - p)
        lines(p, dp, col="purple")
      }
      return(c(sVal = sVal, seSval = seSval))
    }
    
    # Calculate selection coefficients for each locus class
    # Return those separately for each replication of the simulation
    allSval <- NULL
    for (mp in c("main", "purge")){
      for (lc in 1:nClasses){ # Classes 3, 5, and 6 are similar
        allSval <- rbind(allSval, calcSelCoef(lc, mp))
      }
    }
    return(allSval)
  }
}#END calcLocParms

####################################################
# Here is where the analyses start
####################################################
# byIndData organized as:
# 24 times: data.frame with individuals in rows & frequencies for locClasses
# in sets of four columns
byInd <- TRUE
if (byInd){
  outFiles <- list.files(path=paste0(prefix, "AllOut"), pattern="_purge")
  outFiles <- outFiles[-grep("_NoSelf", outFiles)]
  
  print("####### outFiles #######")
  print(paste0(prefix, "AllOut"))
  print(outFiles)
  getGvalFrstPred <- function(bslData){
    getFirstPred <- function(GID){
      bslData$predRec$predict[which(bslData$predRec$predGID == GID)[1]]
    }
    gVal <- bslData$gValue
    basePopID <- bslData$genoRec$basePopID
    popID <- bslData$genoRec$popID
    firstPred <- sapply(1:length(gVal), getFirstPred)
    return(cbind(basePopID, popID, gVal, firstPred))
  }
  calcWgtFreqNdel <- function(byIndDataDF){
    isAdd <- length(grep("Additive", fileName)) > 0
    # The weights look right as of 4 Nov 2019
    if (isAdd) wgt <- c(8, 4, 2, 1, 2) else wgt <- c(8, 4, 2, 1, 2, 2)
    favFreq <- apply(byIndDataDF, 1, function(v24) weighted.mean(v24[1:lenWgt * 4], w=wgt))
    delHomFreq <- apply(byIndDataDF, 1, function(v24) weighted.mean(v24[1:lenWgt * 4 - 3], w=wgt))
    return(cbind(favFreq, delHomFreq))
  }
  
  for (fileName in outFiles){
    # Make memory space
    rm(list=setdiff(ls(), c("calcLocParms", "getGvalFrstPred", "calcWgtFreqNdel", "outFiles", "fileName", "prefix"))); gc()

    outSuff <- substring(fileName, 30)
    
    allOut <- readRDS(paste0(prefix, "AllOut/", fileName))
    byIndData <- lapply(allOut, calcLocParms, analysisToDo="byInd")
    saveRDS(byIndData, file=paste0(prefix, "ByInd/byIndData_", outSuff))
    print("_________________________________________________")
    print(fileName)
    print("_________________________________________________")
    gvfp <- lapply(allOut, getGvalFrstPred)
    ffdhf <- lapply(byIndData, calcWgtFreqNdel)
    byIndParms <- mapply(cbind, gvfp, ffdhf, SIMPLIFY=F)
    
    # Function is poorly named as it calculates the regression coefficient
    # of the favorable allele frequency on the genomic prediction for the
    # main selection event versus on the genotypic value for the purge
    getCorWithFav <- function(bipDF){
      getMainRegCoef <- function(basePop){
        fm <- lm(favFreq ~ firstPred, data=as.data.frame(bipDF[bipDF[,1]==basePop,]))
        return(coef(fm)[2])
      }
      mainCor <- sapply(c(0, 4, 8, 12, 17), getMainRegCoef)
      getPurgeRegCoef <- function(basePop){ # Within-family correlations
        subBipDF <- as.data.frame(bipDF[bipDF[,1]==basePop, c(3,5)])
        ac <- sapply(0:59, function(i) coef(lm(favFreq ~ V1, data=subBipDF[1:10 + i*10,]))[2])
        return(mean(ac))
      }
      purgeCor <- sapply(c(2, 6, 10, 15, 21), getPurgeRegCoef)
      return(c(mainCor, purgeCor))
    }
    allCor <- sapply(byIndParms, getCorWithFav)

    # Similar to the above but it just gives the covariance
    getCovParms <- function(bipDF){
      getMainCov <- function(basePop){
        return(c(var(bipDF[bipDF[,1]==basePop, 5:4]))) # favFreq, firstPred
      }
      mainCov <- sapply(c(0, 4, 8, 12, 17), getMainCov)
      getPurgeCov <- function(basePop){ # Within-family correlations
        subBipDF <- as.data.frame(bipDF[bipDF[,1]==basePop, c(5, 3)]) # favFreq, genoVal
        ac <- sapply(0:59, function(i) c(var(subBipDF[1:10 + i*10, ])))
        return(rowMeans(ac))
      }
      purgeCov <- sapply(c(2, 6, 10, 15, 21), getPurgeCov)
      return(c(mainCov, purgeCov))
    }
    allCov <- sapply(byIndParms, getCovParms)
    
    # Similar to the above but it gives the variance of the independent variable
    getVarIndep <- function(bipDF){
      getMainIndep <- function(basePop){
        return(c(var(bipDF[bipDF[,1]==basePop, 4])))
      }
      mainIndep <- sapply(c(0, 4, 8, 12, 17), getMainIndep)
      getPurgeIndep <- function(basePop){ # Within-family correlations
        subBipDF <- bipDF[bipDF[,1]==basePop, 3]
        ac <- sapply(0:59, function(i) c(var(subBipDF[1:10 + i*10])))
        return(mean(ac))
      }
      purgeIndep <- sapply(c(2, 6, 10, 15, 21), getPurgeIndep)
      return(c(mainIndep, purgeIndep))
    }
    allIndep <- sapply(byIndParms, getVarIndep)
    saveRDS(byIndParms, file=paste0(prefix, "ByInd/byIndAll_", outSuff))
    saveRDS(list(allCor, allCov, allIndep), file=paste0(prefix, "ByInd/byInd_", outSuff))
  }
}#END byInd

# Do this if you just want the class frequencies
# Did not use this in the end
calcClassFreq <- FALSE
if (calcClassFreq){
  locClFreqRep <- lapply(stdSelf, calcLocParms, analysisToDo="locusClass")
  lcf <- locClFreqRep[[1]]
  for (i in 2:24){
    lcf <- lcf + locClFreqRep[[i]]
  }
  View(lcf / 24)
}

###############################################################################
# These were the first analyses
# Overview:
# Take the allele dosage matrix from the simulation and calculate
# freqDelet, freqAllelSt, freqAllelEn, freqLost, freqFixed, freqChng
# for main versus purge, for the different locus effect sizes, and for the
# cycles. And save that.  So they are 2 x (4 or 6) x 5 (maybe 6) matrices
# where the 4 or 6 comes from additive versus dominance.
###############################################################################
# strtEndData organized as:
# 24 times: data.frame for each loc class: vecs for start and end
# freqDelet and freqAllelSt are there to show that there might be a divergence
# between allele frequency of the favorable allele going up but nevertheless
# the frequency of homozygous deleterious recessive genotypes_ also going up
# That would make breeding values increase while genotypic values decrease
calcLocFreq <- FALSE
if (calcLocFreq){
  outFiles <- list.files(path=paste0(prefix, "AllOut"), pattern="_purge")
  remvNoSelf <- grep("_NoSelf", outFiles)
  if (length(remvNoSelf) > 0) outFiles <- outFiles[-remvNoSelf]
  print("####### outFiles #######")
  print(paste0(prefix, "AllOut"))
  print(outFiles)
  for (fileName in outFiles){
    # Make memory space
    rm(list=setdiff(ls(), c("calcLocParms", "outFiles", "fileName", "prefix"))); gc()
    
    stdSelf <- readRDS(paste0(prefix, "AllOut/", fileName))
    strtEndData <- lapply(stdSelf, calcLocParms, analysisToDo="freqDat")
    print("_________________________________________________")
    print(fileName)
    print("_________________________________________________")
    allFreq <- data.frame()
    for (repl in 1:24){
      for (cls in 1:length(strtEndData[[1]])){
        sed <- strtEndData[[repl]][[cls]]
        nr <- nrow(sed)
        colnames(sed)[4:11] <- paste0(rep(c("st", "en"), each=4), colnames(sed)[4:11])
        allFreq <- rbind(allFreq, cbind(data.frame(replication=rep(repl, nr), effSize=rep(locClNames[cls], nr), stringsAsFactors=F), sed))
      }
    }
    
    allFreq$replication <- as.factor(allFreq$replication)
    allFreq$deltaP <- allFreq$enfreq - allFreq$stfreq
    allFreq$pq <- allFreq$stfreq * (1 - allFreq$stfreq)
    allFreq$lost <- ifelse(allFreq$stfreq > 0 & allFreq$enfreq==0, 1, 0)
    allFreq$fixed <- ifelse(allFreq$stfreq < 1 & allFreq$enfreq==1, 1, 0)
    
    freqDelet <- tapply(allFreq$st0, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    freqAllelSt <- tapply(allFreq$stfreq, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    freqAllelEn <- tapply(allFreq$enfreq, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    freqLost <- tapply(allFreq$lost, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    freqFixed <- tapply(allFreq$fixed, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    freqChng <- tapply(allFreq$deltaP, allFreq[, c("replication", "mainPurge", "effSize", "cycle")], mean)
    
    freqDeletSD <- apply(freqDelet, 2:4, sd)
    freqAllelStSD <- apply(freqAllelSt, 2:4, sd)
    freqAllelEnSD <- apply(freqAllelEn, 2:4, sd)
    freqLostSD <- apply(freqLost, 2:4, sd)
    freqFixedSD <- apply(freqFixed, 2:4, sd)
    freqChngSD <- apply(freqChng, 2:4, sd)
    
    freqDelet <- apply(freqDelet, 2:4, mean)
    freqAllelSt <- apply(freqAllelSt, 2:4, mean)
    freqAllelEn <- apply(freqAllelEn, 2:4, mean)
    freqLost <- apply(freqLost, 2:4, mean)
    freqFixed <- apply(freqFixed, 2:4, mean)
    freqChng <- apply(freqChng, 2:4, mean)
    
    doSelCoefProxy <- FALSE
    if (doSelCoefProxy){
      # This loop calculates a proxy for the selection coefficient:
      # The regression of allele frequency change on pq, 
      # where p is the starting allele frequency of the favorable allele.
      # In the standard setting it shows that larger effect alleles have larger
      # selection coefficients
      # Did not use this in the end
      parmsOfInterest <- NULL
      for (mp in c("main", "purge")){
        for (cls in unique(allFreq$effSize)){
          cat(mp, cls, "\n")
          mpClsFreq <- subset(allFreq, allFreq$mainPurge == mp & allFreq$effSize == cls)
          fm <- lmer(deltaP ~ -1 + pq + (-1 + pq | replication), REML=F, data=mpClsFreq)
          sfm <- summary(fm)
          wantedParms <- c(sfm$logLik, sfm$coefficients[-3], attr(sfm$varcor$replication, "stddev"), sfm$sigma)
          # print(summary(fm))
          # print(logLik(fm))
          fm0 <- lm(deltaP ~ -1 + pq, data=mpClsFreq)
          sfm <- summary(fm0)
          wantedParms <- c(wantedParms, logLik(fm0), sfm$coefficients[-(3:4)], sfm$sigma)
          names(wantedParms) <- c("logL_fm", "s_fm", "se_s_fm", "sd_s_rep", "errVar_fm", "logL_fm0", "s_fm0", "se_s_fm0", "errVar_fm")
          # print(summary(fm0))
          # print(logLik(fm0))
          parmsOfInterest <- rbind(parmsOfInterest, wantedParms)
          # plot(mpClsFreq$stfreq, mpClsFreq$deltaP, pch=16, cex=0.8, col=as.integer(substr(mpClsFreq$cycle, start=2, stop=2)), main=paste(mp, cls))
        }
      }
      rownames(parmsOfInterest) <- paste0(rep(c("main", "purge"), each=nClasses), "_", paste0("locCl", 1:nClasses))
    }
    outSuff <- substring(fileName, 30)
    # Suffix St means start, En mean end
    if (doSelCoefProxy){
      saveRDS(list(poi=parmsOfInterest, fDel=freqDelet, fAlSt=freqAllelSt, fAlEn=freqAllelEn, fLost=freqLost, fFixed=freqFixed, fChng=freqChng, fDelSD=freqDeletSD, fAlStSD=freqAllelStSD, fAlEnSD=freqAllelEnSD, fLostSD=freqLostSD, fFixedSD=freqFixedSD, fChngSD=freqChngSD), file=paste0(prefix, "Processed/processed_", outSuff))
    } else{
      saveRDS(list(fDel=freqDelet, fAlSt=freqAllelSt, fAlEn=freqAllelEn, fLost=freqLost, fFixed=freqFixed, fChng=freqChng, fDelSD=freqDeletSD, fAlStSD=freqAllelStSD, fAlEnSD=freqAllelEnSD, fLostSD=freqLostSD, fFixedSD=freqFixedSD, fChngSD=freqChngSD), file=paste0(prefix, "Processed/processed_", outSuff))
    }
    saveRDS(allFreq, file=paste0(prefix, "AllFreq/allFreq_", outSuff))
  }
}#END calcLocFreq

makePlots <- FALSE
if (makePlots){
  #########################################################
  # Process processed
  # Function to plot main versus purge selection gains
  # Works on the first four rows of procVals
  plotProcVals <- function(procVals, plotFileName, parmName){
    plotCol <- sapply(colnames(procVals), function(n) length(grep("Additive", n))*3 + length(grep("MeanSelf", n))*2 + length(grep("_Self", n)))
    plotPch <- sapply(colnames(procVals), function(n) 14 + length(grep("200", n))*3 + length(grep("600", n))*2 + length(grep("1500", n)))
    # Note: the processed is just giving the SD across the 24 reps.  I have NOT
    # divided by sqrt(24).  Now, because sqrt(24) is a constant, it would not have
    # affected any of the weights, so I can just divide by that now.
    xSE <- procVals[3,] / sqrt(24); ySE <- procVals[4,] / sqrt(24)
    xRange <- rbind(procVals[1,] - xSE, procVals[1,] + xSE)
    yRange <- rbind(procVals[2,] - ySE, procVals[2,] + ySE)
    xlim <- range(c(xRange)); ylim <- range(c(yRange))
    xlab <- paste0("Main selection ", parmName, " change")
    ylab <- paste0("Purge selection ", parmName, " change")
    pdf(plotFileName)
    plot(t(procVals[1:2,]), col=plotCol, pch=plotPch, xlab=xlab, ylab=ylab, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlim=xlim, ylim=ylim)
    # Do the stdErr
    for (pt in 1:9){
      lines(xRange[,pt], rep(procVals[2,pt], 2), lwd=2, col=plotCol[pt])
      lines(rep(procVals[1,pt], 2), yRange[,pt], lwd=2, col=plotCol[pt])
    }
    dev.off()
  }
  
  # Inside make plots
  # These processed_ files created above in the section calcLocFreq
  # setwd("/Users/jj332/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/PurgingSimulationScript2019/AnalysisScripts/processedFreq")
  setwd("~/Documents/BioHPC/jj332/PurgeDataScripts/Processed")
  procFiles <- list.files(pattern="processed_")
  procFiles <- c("processed_F1500NORRANAdditiveModelFounder1500SelfAllOut.rds", "processed_F1500NORRANFounderSize1500_MeanSelfAllOut.rds", "processed_F1500NORRANFounderSize1500_SelfAllOut.rds", "processed_F200NORRANAdditiveModelFounder200SelfAllOut.rds", "processed_F200NORRANFounderSize200_MeanSelfAllOut.rds", "processed_F200NORRANFounderSize200_SelfAllOut.rds", "processed_F600NORRANAdditiveModelFounder600SelfAllOut.rds", "processed_F600NORRANStandard_MeanSelfAllOut.rds", "processed_F600NORRANStandard_SelfAllOut.rds")
  parmsToProc <- c("fChng", "fLost")
  # parmsToProc <- c("fChng") # Use this for resPV6
  
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
  dummy <- plotProcVals(procAlFreq, "MainPurgeAlleleFreqChange.pdf", "allele freq.")
  
  resPV6 <- FALSE
  if (do_resPV6){
  # Code to process pV6 which has 1500, 200, and 600 of NoSelf and Self
  pF <- c("processed_F1500NORRANFounderSize1500_NoSelfAllOut.rds", "processed_F1500NORRANFounderSize1500_SelfAllOut.rds", "processed_F200NORRANFounderSize200_NoSelfAllOut.rds", "processed_F200NORRANFounderSize200_SelfAllOut.rds", "processed_F600NORRANStandard_NoSelfAllOut.rds", "processed_F600NORRANStandard_SelfAllOut.rds")
  pV6 <- lapply(pF, procProc)
  tablesTogether <- list(c(1, 3, 5), c(2, 4, 6))
  makeWtMn <- function(tabTog){
    reformTT <- function(tabTog){
      refOneRow <- function(theRow){
        toRet <- NULL
        for(tt in tabTog){
          toRet <- rbind(toRet, pV6[[tt]][[1]][theRow,])
        }
        return(toRet)
      }
      return(lapply(1:4, refOneRow))
    }
    wtMnByCol <- function(theCol){
      wtMn <- function(mp) weighted.mean(refTab[[mp]][, theCol], w=1/refTab[[2+mp]][, theCol])
      sdWtMn <- function(mp) sqrt(1 / sum(1/refTab[[2+mp]][, theCol]))
      return(c(sapply(1:2, wtMn), sapply(1:2, sdWtMn)))
    }
    refTab <- reformTT(tabTog)
    return(sapply(1:5, wtMnByCol))
  }
  resPV6 <- lapply(tablesTogether, makeWtMn)
  saveRDS(resPV6, file="AverageShiftWgtFavAlFreqOverFounderPops.rds")
  tt15_6 <- list(c(1, 5), c(2, 6))
  resPV6_15_6 <- lapply(tt15_6, makeWtMn)
  }
  # These files have genetic gains, i.e. changes in the genotypic values 
  # from one cycle to the next, calculated from the raw simulation outputs
  # Files to be processed here created by "MainVsPurgeGains.R"
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
  dummy <- plotProcVals(procGenVal, "MainPurgeGenValChange.pdf", "genetic value")
  
# In make plots
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
  xlab <- "Genotypic value change"
  ylab <- "Favorable allele freq. change"
  pdf("AlleleFreqVsGenoVal.pdf")
  plot(resGenVal$chng, resAlFr$chng, col=plotCol, pch=plotPch, xlab=xlab, ylab=ylab, cex=1.3, cex.lab=1.3, cex.axis=1.3, xlim=xlim, ylim=ylim)
  # Do the stdErr
  for (pt in 1:18){
    lines(xRange[,pt], rep(resAlFr$chng[pt], 2), lwd=2, col=plotCol[pt])
    lines(rep(resGenVal$chng[pt], 2), yRange[,pt], lwd=2, col=plotCol[pt])
  }
  dev.off()
}
