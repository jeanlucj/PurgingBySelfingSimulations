# WARNING.  These functions assume that the SLiM output is in a variable
# called slimOut, and the program SLiM used is in slimProg.  
# They create tables called mutData and indRaw, which are
# names assumed by downstream functions.  This is bad programming practice
# but it saves moving and duplicating very big files.

# Parse the mutations
parseMut <- function(){
  strtMut <- grep("Mutations", slimOut)[1] + 1
  endMut <- grep("Individuals", slimOut)[1] - 1
  mutRaw <- t(matrix(unlist(strsplit(slimOut[strtMut:endMut], " ")), 9))
  mutData <- data.frame(apply(mutRaw[, -c(1, 3, 7)], 2, as.numeric))
  colnames(mutData) <- c("permID", "basePos", "effect", "degreeDom", "genAppeared", "numGenomes")
  rownames(mutData) <- mutRaw[,1]
  mutData$mutType <- mutRaw[,3]
  mutData$popNum <- mutRaw[,7]
  return(mutData[order(mutData$basePos),])
}

# Parse the individuals
parseInd <- function(){
  strtInd <- grep("Genomes", slimOut)[1] + 1
  nInd <- (length(slimOut) - strtInd + 1) / 2
  indRaw <- strsplit(slimOut[strtInd:length(slimOut)], " ")
  return(lapply(indRaw, function(vec) as.numeric(vec[-(1:2)])))
}

# Calculate own fitness and average fitness of selfed progeny
# For average inbreeding depression, just assume that all loci segregate independently and calculate the mean fitness for 1:2:1.  Log fitness is additive, so I think just calculate the average of 1:2:1 and sum it up with everything else.
# Read the individual
# Create unique list of mutations
# For each mutation calculate 1 + I(Het)*eff*dom + I(Hom)*eff for the individual itself and 1 + I(Het)*(0*0.25 + eff*(0.5*dom + 0.25)) + I(Hom)*eff for its selfed offspring.
# Sum the log of these mutation effects.  Average exponents of these are the fitnesses
calcIndOwnSelfFit <- function(ind){
  indMut <- c(indRaw[[ind*2 - 1]], indRaw[[ind*2]])
  tabMut <- table(indMut)
  eff <- mutData[names(tabMut), "effect"]
  dom <- mutData[names(tabMut), "degreeDom"]
  nMut <- length(eff)
  tabMut <- tabMut[eff != 0]
  dom <- dom[eff != 0]
  eff <- eff[eff != 0]
  nEffMut <- length(eff)
  ownFit <- 1 + cbind(eff*dom, eff)[cbind(1:nEffMut, tabMut)]
  selfFit <- 1 + cbind(eff*(0.25 + 0.5*dom), eff)[cbind(1:nEffMut, tabMut)]
  return(c(sum(log(ownFit)), sum(log(selfFit)), sum(tabMut == 1), nMut, nEffMut))
}

calcOwnSelfFit <- function(indVec){
  fitnesses <- data.frame(t(sapply(indVec, calcIndOwnSelfFit)))
  colnames(fitnesses) <- c("ownFit", "selfedFit", "nEffHet", "nMutations", "nEffMutations")
  return(fitnesses)
}

# Calculate a score that is low for areas of high level of constraint
# It takes into account the rareness of alleles and that rare should be grouped
# which GERP does too, so I call it pseudoGERP
# Assumes the existence of mutData and indRaw in .GlobalEnv
calcPseudoGERP <- function(width=1500, expo=0.5, genomeSize=3e6, hapVec=NULL){
  # Table the size of the genome with zeroes
  # There is a bug that causes names of temp to be in scientific notation
  # The options(scipen=999) is to get rid of that
  makeTable <- function(mutID, genomeSize){
    oldScipen <- options(scipen=999)
    bp <- mutData$basePos[order(as.integer(rownames(mutData)))]
    temp <- table(bp[mutID+1])
    nMutByPos <- integer(genomeSize)
    nMutByPos[as.integer(names(temp))+1] <- temp # One base
    options(scipen=oldScipen)
    return(nMutByPos)
  }
  
  wgtMean <- function(width=1500, expo=0.5, genomeSize=3e6){
    wgtMeanMid <- function(pos) weighted.mean(nMutByPos[pos+posVec]^expo, wgtVec)
    wgtMeanLo <- function(pos){
      pos <- pos + posVec
      weighted.mean(nMutByPos[pos[pos > 0]]^expo, wgtVec[pos > 0])
    }
    wgtMeanHi <- function(pos){
      pos <- pos + posVec
      weighted.mean(nMutByPos[pos[pos <= genomeSize]]^expo, wgtVec[pos <= genomeSize])
    }
    posVec <- 0:width - width/2
    wgtVec <- 1 + width/2 - abs(posVec)
    nLo <- sum(posPoly <= width / 2)
    if (nLo) gerpL <- sapply(posPoly[1:nLo], wgtMeanLo) else gerpL <- NULL
    nHi <- sum(posPoly + width / 2 > genomeSize)
    if (nHi) gerpH <- sapply(posPoly[length(posPoly) - (nHi-1):0], wgtMeanHi) else gerpH <- NULL
    gerpM <- sapply(posPoly[(nLo+1):(length(posPoly) - nHi)], wgtMeanMid)
    data.frame(posPoly=posPoly - 1, gerp=c(gerpL, gerpM, gerpH)) # Zero base
  }

  if (is.null(hapVec)){
    nMutByPos <- makeTable(unlist(indRaw), genomeSize)
  } else{
    nMutByPos <- makeTable(unlist(indRaw[hapVec]), genomeSize)
  }
  
  posPoly <- which(nMutByPos > 0)
  return(wgtMean(width=width, expo=expo, genomeSize=genomeSize))
}

selectMrkOnGERP <- function(gerp, nMrk=480, sdRatio=1, genomeSize=3e6, nIter=10000, plotTrace=F){
  nLoci <- nrow(gerp)
  posPoly <- gerp$posPoly
  gerp <- scale(gerp$gerp)
  upLim <- -0.7
  pNoMove <- 1
  acptStrt <- which(gerp < upLim)
  mrkSet <- sort(sample(acptStrt, nMrk))
  crit <- mean(gerp[mrkSet]) + sdRatio * sd(posPoly[mrkSet[2:nMrk]] - posPoly[mrkSet[1:(nMrk-1)]]) / genomeSize * nMrk
  ac <- crit
  for (i in 1:nIter){
    if (i %% (nIter / 10) == 0){
      upLim <- upLim - 0.1
      acptStrt <- which(gerp < upLim)
      pNoMove <- pNoMove+1
      print(i)
    }
    ms2 <- mrkSet + sample(c(1, rep(2, pNoMove), 3), nMrk, replace=T) - 2
    if (ms2[1] < 1) ms2[1] <- 1
    if (ms2[nMrk] > nLoci) ms2[nMrk] <- nLoci
    ms2[sample(nMrk, 1)] <- sample(acptStrt, 1)
    while(length(unique(ms2)) < nMrk){
      dup <- duplicated(ms2)
      ms2[dup] <- sample(acptStrt, sum(dup))
    }
    ms2 <- sort(ms2)
    crit2 <- mean(gerp[ms2]) + sdRatio * sd(posPoly[ms2[2:nMrk]] - posPoly[ms2[1:(nMrk-1)]]) / genomeSize * nMrk
    if (crit2 < crit){
      mrkSet <- ms2
      crit <- crit2
    }
    ac <- c(ac, crit)
  }
  if (plotTrace) plot(ac)
  return(mrkSet)
}

plotGERPmrk <- function(gerp, limit, mrkSet, genomeSize=3e6, plotMrk=T, plotGenes=T){
  nMrk <- length(mrkSet)
  print(mean(gerp$gerp[mrkSet]))
  print(sd(gerp$posPoly[mrkSet[2:nMrk]] - gerp$posPoly[mrkSet[1:(nMrk-1)]]) / genomeSize * nMrk)
  print(length(unique(mrkSet)))
  gerp <- subset(gerp, posPoly < limit)
  mrkSet <- mrkSet[mrkSet < nrow(gerp)]
  yl <- range(gerp$gerp)
  # Dummy plot to set the coordinates
  plot(gerp$posPoly, gerp$gerp, main=paste("nMrk =", length(mrkSet)), type="n", cex.axis=1.3, xlab="", ylab="")
  if (plotGenes){
    xleft <- 0:(limit %/% 10000 - 1) * 10000 + 8000
    xright <- 1:(limit %/% 10000) * 10000
    rect(xleft, yl[1], xright, yl[2], col="orange", border=NA)
  }
  if (plotMrk){
    for (mrk in mrkSet){
      lines(rep(gerp$posPoly[mrk], 2), yl, col="blue", lwd=2)
    }
  }
  points(gerp$posPoly, gerp$gerp, main=paste("nMrk =", length(mrkSet)), pch=16, cex=0.8)
}

if (FALSE){
  s <- Sys.time()
  test <- selectMrkOnGERP(gerp, nMrk=1200, sdRatio=2, nIter=100000, plotTrace=T)
  print(Sys.time() - s)
  plotGERPmrk(gerp, 60000, test)
}

# Make a BSL environment out of the SLiM output
# Make the founder haplotypes
# hapVec is a vector of the haplotypes you want
makeFounderHaps <- function(hapVec, mutData){
  genoVecFromInd <- function(mutVec, maxPos){
    genoVec <- integer(maxPos)
    genoVec[rownames(mutData) %in% mutVec] <- 1
    return(genoVec)
  }
  return(t(sapply(indRaw[hapVec], genoVecFromInd, maxPos=nrow(mutData))))
}

# Use SLiM output to create an environment like what comes from defineSpecies
makeBSLspecies <- function(hapVec, nMrk, nQTL, oracle, selMrk=F, allQ=T){
  # Make the map data
  # Chr and Pos are based on known parameters given to SLiM
  # Figure out chromosome size
  recLine <- slimProg[grep("ends", slimProg)[1]]
  recVec <- unlist(strsplit(recLine, "[(,) ]"))
  recVec <- suppressWarnings(as.integer(recVec))
  nBasePairsPerChr <- max(recVec, na.rm=T) + 1
  # Figure out recombination rate
  rateLine <- slimProg[grep("rates", slimProg)[1]]
  rateVec <- unlist(strsplit(rateLine, "[(,) ]"))
  rateVec <- suppressWarnings(as.numeric(rateVec))
  cMperBasePair <- mean(rateVec, na.rm=T) * 100
  # Figure out which mutations are neutral and which causal
  genElemLines <- slimProg[grep("initializeGenomicElementType", slimProg)]
  neutMutNames <- NULL
  # Neutral mutations are always the first in a genetic element
  # Function returns the causal mut names but assigns neut mut names to that vector
  getNeutMutNames <- function(genElemLine){
    mutNames <- unlist(strsplit(genElemLine, "[(,) ]"))
    mutNames <- mutNames[substr(mutNames, 1, 1) == "m"]
    neutMutNames <<- c(neutMutNames, mutNames[1])
    return(mutNames[-1])
  }
  causMutNames <- unlist(sapply(genElemLines, getNeutMutNames))
  
  # Subset mutData so only mutations in the haplotypes of hapVec
  # Should be polymorphic in the final dataset
  subMut <- unique(unlist(indRaw[hapVec]))
  mutData <- mutData[as.character(subMut),] # indRaw uses mutID not idx
  mutData <- mutData[order(mutData$basePos),]
  # Sample the QTL to be at all sites 
  # WARNING 300 QTL hard coded here
  # WARNING particular SLiM program hard coded here
  effectivePos <- which(mutData$mutType %in% causMutNames)
  effPosBP <- mutData$basePos[effectivePos]
  markerPos <- which(mutData$mutType %in% neutMutNames)
  if (allQ){
    ep <- NULL
    for (q in 0:299){
      effCan <- effectivePos[9e3 <= effPosBP - q*1e4 & effPosBP - q*1e4 < 1e4]
      if (length(effCan) > 0) ep <- c(ep, sample(effCan, 1))
    }
    effectivePos <- ep
  } else{
    effectivePos <- sample(effectivePos, nQTL)
  }
  if (oracle){
    markerPos <- c(effectivePos, sample(markerPos, nMrk - length(effectivePos)))
  } else{
    if (selMrk){ # select markers or QTL on pseudoGERP scores
      # NOTE selectMrkOnGERP uses a data.frame index by position
      # NOT by mutID
      markerPos <- selMrkGERP
      gerp <- gerp[gerp$posPoly %in% mutData$basePos,]
      markerPos <- gerp$posPoly[markerPos]
      markerPosDup <- markerPos[duplicated(markerPos)]
      markerPos <- markerPos[!duplicated(markerPos)]
      # Retranslate back onto the mutData index
      basePos <- mutData$basePos
      basePosDupIdx <- which(duplicated(basePos))
      basePosDup <- basePos[duplicated(basePos)]
      basePosIdx <- which(!duplicated(basePos))
      basePos <- basePos[!duplicated(basePos)]
      markerPos <- c(basePosIdx[basePos %in% markerPos], basePosDupIdx[basePosDup %in% markerPosDup])
    } else{ # Randomly pick markers that are neutral
      markerPos <- sample(markerPos, nMrk)
    }
  }
  loci <- sort(union(markerPos, effectivePos))
  mutData <- mutData[loci,]
  markerPos <- which(loci %in% markerPos)
  effectivePos <- which(loci %in% effectivePos)
  
  Chr <- mutData$basePos %/% nBasePairsPerChr + 1
  Pos <- mutData$basePos %% nBasePairsPerChr * cMperBasePair
  map <- data.frame(Chr=Chr, Pos=Pos)
  
  effectID <- 1:length(effectivePos)
  effects <- matrix(mutData[effectivePos, "effect"], ncol=1)
  actionType <- mutData[effectivePos, "degreeDom"]
  mapData <- list(map=map, markerPos=markerPos, effectID=effectID, effectivePos=effectivePos, effects=effects, actionType=actionType, domModel="Partial")
  
  founderHaps <- makeFounderHaps(hapVec, mutData)
  return(list(mapData=mapData, founderHaps=founderHaps))
}#END makeBSLspecies

makeSimEnv <- function(hapVec, nMrk, nQTL, nCore=1, oracle=F, selMrk=F){
  sims <- lapply(hapVec, makeBSLspecies, nMrk, nQTL, oracle, selMrk)
  nSim <- length(hapVec)
  toRemove <- c(setdiff(ls(), c("sims", "nSim", "nCore")), "toRemove")
  rm(list=toRemove)
  return(environment())
}

# All the chromosomes in addEnv will be added into simEnv$sims
addSimEnvChr <- function(simEnv, addEnv){
  for (sim in 1:min(length(simEnv$sims), length(addEnv$sims))){ # Hopefully, they are the same length...
    nLoc <- nrow(simEnv$sims[[sim]]$mapData$map)
    
    addMap <- addEnv$sims[[sim]]$mapData$map
    addMap$Chr <- addMap$Chr + max(simEnv$sims[[sim]]$mapData$map$Chr)
    simEnv$sims[[sim]]$mapData$map <- rbind(simEnv$sims[[sim]]$mapData$map, addMap)
    
    simEnv$sims[[sim]]$mapData$markerPos <- c(simEnv$sims[[sim]]$mapData$markerPos, addEnv$sims[[sim]]$mapData$markerPos + nLoc)
    simEnv$sims[[sim]]$mapData$effectivePos <- c(simEnv$sims[[sim]]$mapData$effectivePos, addEnv$sims[[sim]]$mapData$effectivePos + nLoc)
    
    simEnv$sims[[sim]]$mapData$effectID <- c(simEnv$sims[[sim]]$mapData$effectID, addEnv$sims[[sim]]$mapData$effectID + max(simEnv$sims[[sim]]$mapData$effectID))
    
    simEnv$sims[[sim]]$mapData$effects <- rbind(simEnv$sims[[sim]]$mapData$effects, addEnv$sims[[sim]]$mapData$effects)
    simEnv$sims[[sim]]$mapData$actionType <- c(simEnv$sims[[sim]]$mapData$actionType, addEnv$sims[[sim]]$mapData$actionType)
    
    simEnv$sims[[sim]]$founderHaps <- cbind(simEnv$sims[[sim]]$founderHaps, addEnv$sims[[sim]]$founderHaps)
  }
  return(max(addMap$Chr))
}

# Combine simulations into on environment so they can be processed in parallel
combineSimEnvSims <- function(simEnv, addEnv){
  simEnv$sims <- c(simEnv$sims, addEnv$sims)
  simEnv$nSim <- simEnv$nSim + addEnv$nSim
  return(simEnv$nSim)
}

# Watterson's estimator of Ne
# This one needs a lot of work
calcWattersonsNe <- function(){
  sampSize <- 100
  correction <- sum(1 / 1:sampSize)
  neutMut <- as.numeric(rownames(mutData)[mutData$mutType == "m1"])
  nEstimates <- 500
  neEst <- NULL
  for (est in 1:nEstimates){
    samp <- sample(length(indRaw), sampSize)
    mutInSamp <- intersect(unique(unlist(indRaw[samp])), neutMut)
    neEst <- c(neEst, length(mutInSamp) / correction / 4 / (6e6 * (0.85*(1e-7) + 0.15*(2.5e-8))))
  }
}
