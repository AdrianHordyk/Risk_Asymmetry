


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

betaVec <- exp(seq(-0.5, 0.5, length.out=9))

# Run historical 
message("Calculate depletion")
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
OM@D <- rep(dep, 2)

runOM <- OM
nsim <- runOM@nsim

temp <- list()
for (x in 1:length(betaVec)) {
  beta <- betaVec[x]
  message(x, "/", length(betaVec),  " -  beta = ", round(beta,2))
  
  if (length(names(OM@cpars)) > 0) {
    cparnsim <- dim(OM@cpars[[1]])[1]
    runOM@cpars$beta <- rep(beta, cparnsim)
  } else {
    runOM@cpars$beta <- rep(beta, nsim)
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
  outDF$var <- beta
  outDF$Name <- MSE@Name
  temp[[x]] <- outDF
  
}

DF <- do.call('rbind', temp)
Name <- gsub(" ", "_", OM@Name)
flname <- paste0("betabias_", Name, "_", Dep, ".rdata")
saveRDS(DF, file.path(Resultsdir, flname))




# betaVec <- seq(0.5, 2, length.out = 7)

# 
# for (Dep in DepVec) {
#   for (beta in betaVec) {
#     runOM <- OM 
#     Hist <- runMSE(runOM, Hist=TRUE)
#     dep <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * Dep 
#     runOM@D <- rep(dep, 2)
#     if (length(names(OM@cpars)) > 0) {
#       cparnsim <- dim(OM@cpars[[1]])[1]
#       runOM@cpars$beta <- rep(beta, cparnsim)
#     } else {
#       runOM@cpars$beta <- rep(beta, nsim)
#     }
#     
#     message(OM@Name, " - ", Dep)
#     MSE <- runMSE(runOM, MPs='SCA_MP', parallel = TRUE)
#     
#     omdat <- MSE@OM
#     omdat$sim <- 1:MSE@nsim
#     obsdat <-MSE@Obs
#     obsdat$sim <- 1:MSE@nsim
#     
#     output <- data.frame(B_BMSY=as.vector(unlist(MSE@B_BMSY[,1,])),
#                          F_FMSY=as.vector(unlist(MSE@F_FMSY[,1,])),
#                          Catch=as.vector(unlist(MSE@C[,1,])),
#                          TAC=as.vector(unlist(MSE@TAC[,1,])),
#                          sim=1:MSE@nsim, Years=rep(1:MSE@proyears, each=MSE@nsim))
#     
#     outDF <- left_join(output, omdat, by='sim')
#     outDF <- left_join(outDF, obsdat, by='sim')
#     outDF$dep <- dep
#     outDF$var <- beta
#     outDF$Name <- MSE@Name
#     
#     Name <- gsub(" ", "_", OM@Name)
#     flname <- paste0("Beta_", Name, "_", beta, "_", Dep, ".rdata")
#     saveRDS(outDF, file.path(Resultsdir, flname))
#   }
# }
