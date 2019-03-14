

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

A50 <- OM@SizeLimSD[1]
ageMVec <-(1+ seq(-0.5, 0.5, length.out = 9))

assumeA50s <- A50 * ageMVec

# Calculate bias in age-of-maturity
Ages <- 1:OM@maxage
trueL50 <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(A50-mean(OM@t0))))
assumeL50s <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(assumeA50s-mean(OM@t0))))

L50biasvec <- assumeL50s/trueL50


# Run historical 
message("Calculate depletion")
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
OM@D <- rep(dep, 2)

runOM <- OM
nsim <- runOM@nsim

temp <- list()
for (x in 1:length(L50biasvec)) {
  L50bias <- L50biasvec[x]
  message(x, "/", length(L50biasvec),  " -  L50bias = ", round(L50bias,2))
  
  if (length(names(OM@cpars)) > 0) {
    cparnsim <- dim(OM@cpars[[1]])[1]
    runOM@cpars$lenMbias <- rep(L50bias, cparnsim)
  } else {
    runOM@cpars$lenMbias <- rep(L50bias, nsim)
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
  outDF$var <- ageMVec[x]
  outDF$Name <- MSE@Name
  temp[[x]] <- outDF
  
}

DF <- do.call('rbind', temp)
Name <- gsub(" ", "_", OM@Name)
flname <- paste0("AgeMbias_", Name, "_", Dep, ".rdata")
saveRDS(DF, file.path(Resultsdir, flname))


# AgeMVec <- seq(-2, 2, length.out = 5)
# 
# # True AgeM 
# BaseMage <- OM@cpars$Mat_age[1,,1]
# baseAm <- min(which(BaseMage >=0.5))
# 
# ageMVec <-(1+ seq(-0.5, 0.5, length.out = 9))
# 
# OM_base <- OM 


# OM_Low <- OM # low ageM
# lowAgeM <- c(BaseMage, 1, 1)
# lowAgeM <- matrix(lowAgeM[3:length(lowAgeM)], nrow=OM@nsim, ncol=OM@maxage, byrow=TRUE)
# OM_Low@cpars$Mat_age <- array(lowAgeM, dim=dim(OM@cpars$Mat_age))
# 
# 
# OM_High <- OM # high ageM 
# highAgeM <- c(0,0, BaseMage)
# highAgeM <- matrix(highAgeM[1:(length(highAgeM)-2)], nrow=OM@nsim, ncol=OM@maxage, byrow=TRUE)
# OM_High@cpars$Mat_age <- array(highAgeM, dim=dim(OM@cpars$Mat_age))



# #### DEBUG ####
# useOM <- OM
# A50 <- min(which(useOM@cpars$Mat_age[1,,1] >=0.5))
# assumeA50s <- AgeMVec + A50
# 
# Ages <- 1:OM@maxage
# trueL50 <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(A50-mean(OM@t0))))
# assumeL50s <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(assumeA50s-mean(OM@t0))))
# 
# L50biasvec <- assumeL50s/trueL50
# 
# useOM@cpars$lenMbias <- rep(L50biasvec[1], useOM@nsim)
# 
# M1 <- runMSE(useOM, Hist=TRUE)
# 
# iVB(mean(useOM@t0), mean(useOM@K), mean(useOM@Linf), M1$SampPars$L50[1])
# iVB(mean(useOM@t0), mean(useOM@K), mean(useOM@Linf), M1$Data@L50[1])
# 
# M1$Data@L50/M1$SampPars$L50
# 
# file <- tempfile()
# saveRDS(useOM, file)
# 
# OM <- readRDS("C:\\Users\\Adrian\\AppData\\Local\\Temp\\RtmpU3ivZ9\\file1bc855881c42")
# ### End Debug ####

# 
# runAm <- function(OM, Dep, L50bias, parallel=TRUE) {
#   runOM <- OM 
#   
#   Hist <- runMSE(runOM, Hist=TRUE)
#   dep <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * Dep 
#   
#   message(OM@Name)
#   if (length(names(OM@cpars)) > 0) {
#     cparnsim <- dim(OM@cpars[[1]])[1]
#     runOM@cpars$lenMbias <- rep(L50bias, cparnsim)
#   } else {
#     runOM@cpars$lenMbias <- rep(L50bias, nsim)
#   }
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
#   outDF$var <- L50bias
#   outDF$Name <- MSE@Name
#   outDF
# }
# 
# 
# 
# A50 <- OM@SizeLimSD[1]
# assumeA50s <- ageMVec * A50
# 
# Test <- "AgeMbias"
# # Base case 
# message("Base Case")
# useOM <- OM_base
# Case <- "Base"
# 
# # Calculate bias in age-of-maturity
# Ages <- 1:OM@maxage
# trueL50 <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(A50-mean(OM@t0))))
# assumeL50s <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(assumeA50s-mean(OM@t0))))
# 
# L50biasvec <- assumeL50s/trueL50
# 
# for (Dep in DepVec) {
#   temp <- lapply(1:length(L50biasvec), function(x) runAm(useOM, Dep, L50biasvec[x]))
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
# }
# 

# 
# # Low age-M 
# useOM <- OM_Low
# Case <- "Low"
# 
# # Calculate bias in age-of-maturity
# A50_low <- useOM@SizeLimSD[1] - 2
# assumeA50s <- A50_low + c(0, 2, 4)
# Ages <- 1:OM@maxage
# trueL50 <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(A50-mean(OM@t0))))
# assumeL50s <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(assumeA50s-mean(OM@t0))))
# 
# L50biasvec <- assumeL50s/trueL50
# L50biasvec <- c(min(L50biasvec), median(L50biasvec), max(L50biasvec))
# 
# for (Dep in DepVec) {
#   temp <- lapply(1:length(L50biasvec), function(x) runAm(useOM, Dep, L50biasvec[x]))
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
#   
# }
# 
# 
# # High age-M 
# useOM <- OM_High
# Case <- "High"
# # Calculate bias in age-of-maturity
# A50_high <- useOM@SizeLimSD[1] + 2
# assumeA50s <- A50_high + c(-4, -2, 0)
# Ages <- 1:OM@maxage
# trueL50 <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(A50-mean(OM@t0))))
# assumeL50s <- mean(OM@Linf) * (1-exp(-mean(OM@K)*(assumeA50s-mean(OM@t0))))
# 
# L50biasvec <- assumeL50s/trueL50
# L50biasvec <- c(min(L50biasvec), median(L50biasvec), max(L50biasvec))
# 
# for (Dep in DepVec) {
#   temp <- lapply(1:length(L50biasvec), function(x) runAm(useOM, Dep, L50biasvec[x]))
#   DF <- do.call('rbind', temp)
#   Name <- gsub(" ", "_", OM@Name)
#   flname <- paste0(Name, "_", Case, "_", Test, "_", Dep, ".rdata")
#   saveRDS(DF, file.path(Resultsdir, flname))
#   
# }

