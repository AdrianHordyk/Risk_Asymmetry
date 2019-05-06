
args <- commandArgs(trailingOnly = TRUE)
if (length(args)>0) {
  omnum <- as.numeric(args[1])  
  Dep <- as.numeric(args[2])  
} else {
  omnum <- 1
  Dep <- 0.5
}

source("R/control.R")

OM <- MakeOM(omnum, nsim) # load OM 

MbiasVec <- exp(seq(-0.5, 0.5, length.out = 9)) 

# Run historical 
message("Calculate depletion")
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
OM@D <- rep(dep, 2)

runOM <- OM
nsim <- runOM@nsim

temp <- list()
for (x in 1:length(MbiasVec)) {
  Mbias <- MbiasVec[x]
  message(x, "/", length(MbiasVec),  " -  Mbias = ", round(Mbias,2))
  
  if (length(names(OM@cpars)) > 0) {
    cparnsim <- dim(OM@cpars[[1]])[1]
    runOM@cpars$Mbias <- rep(Mbias, cparnsim)
  } else {
    runOM@cpars$Mbias <- rep(Mbias, nsim)
  }

  MSE <- runMSE(runOM, MPs='SCA_MP', parallel = TRUE)
  
  omdat <- MSE@OM
  omdat$sim <- 1:MSE@nsim
  obsdat <-MSE@Obs
  obsdat$sim <- 1:MSE@nsim
  
  output <- data.frame(B_BMSY=as.vector(unlist(MSE@B_BMSY[,1,])),
                       F_FMSY=as.vector(unlist(MSE@F_FMSY[,1,])),
                       Catch=as.vector(unlist(MSE@C[,1,])),
                       TAC=as.vector(unlist(MSE@TAC[,1,])),
                       sim=1:MSE@nsim, Years=rep(1:MSE@proyears, each=MSE@nsim))
  
  outDF <- left_join(output, omdat, by='sim')
  outDF <- left_join(outDF, obsdat, by='sim')
  outDF$dep <- dep
  outDF$var <- Mbias
  outDF$Name <- MSE@Name
  temp[[x]] <- outDF
  
}

DF <- do.call('rbind', temp)
Name <- gsub(" ", "_", OM@Name)
flname <- paste0("Mbias_", Name, "_", Dep, ".rdata")
saveRDS(DF, file.path(Resultsdir, flname))

  

