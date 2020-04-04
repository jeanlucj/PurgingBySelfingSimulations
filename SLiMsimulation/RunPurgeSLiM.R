setwd("~/Google Drive/Teaching/People/MohamedSomoIbrahim/PurgeSimulations/SLiM")

for (slimSim in 1:24){
  simProg <- readLines("PurgeSimRedo20191031.slim")
  saveLine <- grep("#", simProg, fixed=T)
  simProg[saveLine] <- sub("#", slimSim, simProg[saveLine])
  writeLines(simProg, "PurgeSimRedo20191031op.slim")
  system(paste0("slim PurgeSimRedo20191031op.slim > PurgeSimOut", slimSim, ".txt"), intern=F, wait=T)
  print(paste0("slimSim ", slimSim))
}