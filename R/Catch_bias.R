
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


# Run historical 
message("Calculate depletion")
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
OM@D <- rep(dep, 2)

runOM <- OM
nsim <- runOM@nsim

CbiasVec <- 1 + seq(-0.5, to=0.5, length.out=9)

nyears <- OM@nyears
pyears <- OM@proyears
Cyears <- 1:ceiling(2/3*nyears)

temp <- list()
for (x in 1:length(CbiasVec)) {
  Cbias <- CbiasVec[x]
  message(x, "/", length(CbiasVec),  " -  Cbias = ", round(Cbias,2))
  
  MSE <- runMSE(runOM, MPs='SCA_MP', parallel = TRUE, control=list(yrs=Cyears, Cbias_yr=Cbias))
  
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
  outDF$var <- Cbias
  outDF$Name <- MSE@Name
  temp[[x]] <- outDF
  
}

DF <- do.call('rbind', temp)
Name <- gsub(" ", "_", OM@Name)
flname <- paste0("Cbias_", Name, "_", Dep, ".rdata")
saveRDS(DF, file.path(Resultsdir, flname))

# 
# # Run historical 
# message("Calculate depletion")
# Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
# dep <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * Dep
# OM@D <- rep(dep, 2)
# 
# runOM <- OM
# nsim <- runOM@nsim
# 
# Hist <- runMSE(OM, Hist=TRUE)
# Deps <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * DepVec 
# 
# 
# 
# runCatch <- function(OM, dep, Cbias, Cyears, parallel=TRUE) {
#   runOM <- OM 
#   message(OM@Name)
# 
#   runOM@D <- rep(dep, 2)
#   
#   MSE <- runMSE(runOM, MPs='SCA_MP', parallel = parallel, control=list(yrs=Cyears, Cbias_yr=Cbias))
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
#   outDF$var <- Cbias
#   outDF$Name <- MSE@Name
#   outDF
# }
# 
# 
# 
# # --- Catch bias - historical ----
# nyears <- OM@nyears
# pyears <- OM@proyears
# yrs <- 1:ceiling(2/3*nyears)
# 
# count <- 0; templist <- list()
# for (dep in Deps) {
#   for (Cbias in CbiasVec) {
#     count <- count + 1
#     templist[[count]] <- runCatch(OM, dep, Cbias, Cyears=yrs, parallel=TRUE)
#   }
# }
# Results <- do.call('rbind', templist)
# Name <- gsub(" ", "_", OM@Name)
# flname <- paste0(Name, "_Cbias_Hist.rdata")
# saveRDS(Results, file.path(Resultsdir, flname))
#  


 
# # --- Catch bias - recent ----
# 
# yrs <- 1:ceiling(1/3*nyears)
# 
# count <- 0; templist <- list()
# for (dep in Deps) {
#   for (Cbias in CbiasVec) {
#     count <- count + 1
#     templist[[count]] <- runCatch(OM, dep, Cbias, Cyears=yrs, parallel=TRUE)
#   }
# }
# Results <- do.call('rbind', templist)
# Name <- gsub(" ", "_", OM@Name)
# flname <- paste0(Name, "_Cbias_Recent.rdata")
# saveRDS(Results, file.path(Resultsdir, flname))
# 
# # --- Catch bias - all ----
# 
# yrs <- 1:(nyears+pyears)
# 
# count <- 0; templist <- list()
# for (dep in Deps) {
#   for (Cbias in CbiasVec) {
#     count <- count + 1
#     templist[[count]] <- runCatch(OM, dep, Cbias, Cyears=yrs, parallel=TRUE)
#   }
# }
# Results <- do.call('rbind', templist)
# Name <- gsub(" ", "_", OM@Name)
# flname <- paste0(Name, "_Cbias_All.rdata")
# saveRDS(Results, file.path(Resultsdir, flname))
# 
