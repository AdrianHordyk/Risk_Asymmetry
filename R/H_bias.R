
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

CR <- h2CR(OM@h[1])

biasVec <- seq(-0.5, 0.5, length.out = 9)
HbiasVec <-  CR2h(CR * (1+biasVec)) - OM@h[1] 


# Run historical 
message("Calculate depletion")
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
OM@D <- rep(dep, 2)

runOM <- OM
nsim <- runOM@nsim

temp <- list()
for (x in 1:length(HbiasVec)) {
  Hbias <- HbiasVec[x]
  message(x, "/", length(HbiasVec),  " -  Hbias = ", round(Hbias,2))
  OMHval <- mean(OM@h)
  if (length(names(OM@cpars)) > 0) {
    cparnsim <- dim(OM@cpars[[1]])[1]
    runOM@cpars$hsim <- rep(Hbias+OMHval, cparnsim)
  } else {
    runOM@cpars$hsim <- rep(Hbias+OMHval, nsim)
  }
  runOM@cpars$hsim[runOM@cpars$hsim>0.99] <- 0.99
  
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
  outDF$var <- Hbias
  outDF$Name <- MSE@Name
  temp[[x]] <- outDF
  
}

DF <- do.call('rbind', temp)
Name <- gsub(" ", "_", OM@Name)
flname <- paste0("Hbias_", Name, "_", Dep, ".rdata")
saveRDS(DF, file.path(Resultsdir, flname))

# 
# 
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args)>0) {
#   omnum <- as.numeric(args[1])  
# } else {
#   omnum <- 2
# }
# 
# source("R/control.R")
# 
# # Calculate initial depletion as fraction of BMSY
# OM <- readRDS(file.path(OMdir, OMFiles[omnum]))
# nsim <- OM@nsim
# 
# # OM@maxF <- 1.5 
# 
# 
# 
# 
# # True Hs 
# baseH <- mean(OM@h)
# OM_base <- OM 
# # OM_Low <- OM # low h
# # OM_Low@h <- rep(min(HbiasVec) + baseH,2)
# # 
# # OM_High <- OM # high h
# # OM_High@h <- rep(max(HbiasVec) + baseH,2)
# 
# # assumedhs <- mean(OM@h) + c(HbiasVec[1], 0, HbiasVec[length(HbiasVec)])
# Test <- "hbias"
# 
# runH <- function(OM, Dep, Hbias, parallel=TRUE) {
#   runOM <- OM 
#   runOM@h[runOM@h>=1] <- 0.99
#   Hist <- runMSE(runOM, Hist=TRUE)
#   dep <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * Dep 
#   
#   message(OM@Name)
#   OMHval <- mean(OM@h)
#   if (length(names(OM@cpars)) > 0) {
#     cparnsim <- dim(OM@cpars[[1]])[1]
#     runOM@cpars$hsim <- rep(Hbias+OMHval, cparnsim)
#   } else {
#     runOM@cpars$hsim <- rep(Hbias+OMHval, nsim)
#   }
#   
#   runOM@cpars$hsim[runOM@cpars$hsim>0.99] <- 0.99
#   
#   runOM@D <- rep(dep, 2)
#   
#   MSE <- runMSE(runOM, MPs='SCA_MP', parallel = parallel)
#   
#   omdat <- MSE@OM
#   omdat$sim <- 1:MSE@nsim
#   obsdat <-MSE@Obs
#   obsdat$sim <- 1:MSE@nsim
#   
#   output <- data.frame(B_BMSY=as.vector(unlist(MSE@B_BMSY[,1,])),
#                        F_FMSY=as.vector(unlist(MSE@F_FMSY[,1,])),
#                        Catch=as.vector(unlist(MSE@C[,1,])),
#                        TAC=as.vector(unlist(MSE@TAC[,1,])),
#                        sim=1:MSE@nsim, Years=rep(1:MSE@proyears, each=MSE@nsim))
#   
#   outDF <- left_join(output, omdat, by='sim')
#   outDF <- left_join(outDF, obsdat, by='sim')
#   outDF$dep <- dep
#   outDF$var <- Hbias
#   outDF$Name <- MSE@Name
#   outDF
# }
# 
# # Base case
# message("Base Case")
# useOM <- OM_base
# Test <- "hbias"
# Case <- "Base"
# for (Dep in DepVec) {
#   message(OM@Name)
#   message(Dep)
#   temp <- list()
#   for (x in 1:length(HbiasVec)) {
#     message(HbiasVec[x])
#     temp[[x]] <-runH(useOM, Dep, HbiasVec[x])
#   }
#   
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
# }
# 


# # Low h
# message("Low Case")
# useOM <- OM_Low
# 
# Case <- "Low"
# 
# HbiasVec2 <- assumedhs - mean(useOM@h)
# 
# for (Dep in DepVec) {
#   temp <- lapply(1:length(HbiasVec2), function(x) runH(useOM, Dep, HbiasVec2[x]))
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
# 
# }
# 
# # # High h 
# message("High Case")
# useOM <- OM_High
# Case <- "High"
# HbiasVec2 <-  assumedhs - mean(useOM@h)
# 
# for (Dep in DepVec) {
#   temp <- lapply(1:length(HbiasVec2), function(x) runH(useOM, Dep, HbiasVec2[x]))
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
# 
# }
