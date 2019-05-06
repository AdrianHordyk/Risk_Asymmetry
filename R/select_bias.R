
devtools::install_github('tcarruth/MSEtool', ref = "fix_dome")


source("R/control.R")

Vmaxage <- 0.65
AssumedVmax <- Vmaxage * (1 + seq(-0.5, 0.5, length.out=9))

for (omnum in 1:4) {
  message("OM:", omnum)
  OM <- MakeOM(omnum, nsim) # load OM 
  om <- OM 
  om@nsim <- 2
  Hist <- runMSE(om, Hist=TRUE)
  afs <- min(which(Hist@AtAge$Select[1,,1] >=0.99))
  as5 <- min(which(Hist@AtAge$Select[1,,1] >=0.05))
  
  srs <- (om@maxage - afs) / ((-log(Vmaxage,2))^0.5) # selectivity parameters are constant for all years
  srs[!is.finite(srs)] <- Inf
  sls <- (afs - as5) /((-log(0.05,2))^0.5)
  selage <- DLMtool:::getsel(lens=1:om@maxage, lfs=afs, sls=sls, srs=srs)
  
  V <- array(selage, dim=c(om@maxage, om@nsim, om@nyears+om@proyears)) 
  V <- aperm(V, c(2,1,3))
  OM@cpars$V <- V
  
  # Create MPs 
  for (x in seq_along(AssumedVmax)) {
    
    name <- paste('SCA_MP', AssumedVmax[x], sep="_")
    if (AssumedVmax[x]<1) {
      temp <- make_MP(SCA, HCR_MSY, fix_h=TRUE, fix_tau=TRUE, 
                      vulnerability = "dome", Vmaxage=AssumedVmax[x])
    } else {
      temp <- make_MP(SCA, HCR_MSY, fix_h=TRUE, fix_tau=TRUE, 
                      vulnerability = "logistic")
    }
    assign(name, temp)
  }
  
  mps <- avail('MP')
  testMPs <- mps[grepl('SCA_MP_', mps)]
  
  # Run historical 
  message("Calculate depletion")
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  dep <- round(Hist@Ref$SSBMSY_SSB0[1],3) * Dep
  OM@D <- rep(dep, 2)
  
  runOM <- OM
  nsim <- runOM@nsim
  
  MSE <- runMSE(runOM, MPs=testMPs, parallel = TRUE)
  
  omdat <- MSE@OM
  omdat$sim <- 1:MSE@nsim
  obsdat <-MSE@Obs
  obsdat$sim <- 1:MSE@nsim
  
  output <- data.frame(B_BMSY=as.vector(unlist(MSE@B_BMSY)),
                       F_FMSY=as.vector(unlist(MSE@F_FMSY)),
                       Catch=as.vector(unlist(MSE@C)),
                       TAC=as.vector(unlist(MSE@TAC)),
                       sim=1:MSE@nsim, Years=rep(1:MSE@proyears, each=MSE@nMPs*MSE@nsim),
                       MP=rep(MSE@MPs, each=MSE@nsim))
  
  outDF <- left_join(output, omdat, by='sim')
  outDF <- left_join(outDF, obsdat, by='sim')
  outDF$dep <- dep
  outDF$Name <- MSE@Name
  
  DF <- outDF
  Name <- gsub(" ", "_", OM@Name)
  flname <- paste0("select_", Name, "_", Dep, ".rdata")
  saveRDS(DF, file.path(Resultsdir, flname))
}











# 
# 
# 
# 
# Create Figure
# SelectList <- lapply(1:length(AssumedVmax), function (x)
#   selA(as5,afs, Vmaxage=AssumedVmax[x], OM@maxage, OM@nsim, OM))
# 
# OMlist <- list()
# TrueVmax <- Vmaxage
# for (x in seq_along(TrueVmax)) {
#   temp <- OM
#   temp@cpars$V <- SelectList[[x]]
#   name <- paste('OM', TrueVmax[x], sep="_")
#   OMlist[[name]] <- temp
# }
# 
# # Create a plot of selectivity - or write out data.frame
# library(purrr)
# List <- purrr::map(SelectList, function(x) x[1,,1])
# df <- do.call("rbind", List) %>% t() %>% data.frame()
# colnames(df) <- AssumedVmax
# df <- tidyr::gather(df)
# df$Age <- 1:OM@maxage
# df$bias <- as.numeric(df$key)/0.65 -1
# 
# text.size <- 5
# axis.size <- 4
# leg.size <- 5
# 
# df$abs_bias <- abs(df$bias)
# df$bias_direction <- "Negative"
# df$bias_direction[df$bias>0] <- 'Positive'
# df$bias_direction[df$bias==0] <- 'Unbiased'
# df$bias_direction <- factor(df$bias_direction, levels=c('Unbiased', 'Positive', 'Negative'), ordered=TRUE)
# 
# pout <- ggplot(df, aes(Age, value, color=bias_direction, group=bias, linetype=as.factor(abs_bias))) + 
#   geom_line() +
#   scale_color_discrete() +
#   theme_classic() +
#   labs(linetype='Absolute bias Vmax', x="Age (year)", y="Selectivity", color="Bias") +
#   theme(axis.title = element_text(size=text.size),
#         axis.text = element_text(size=axis.size),
#         legend.text = element_text(size=leg.size),
#         legend.key.size = unit(0.6, 'lines'),
#         legend.title = element_text(size=text.size))
# 
# ggsave("Figures/Figure2.png", dpi=600, width=80, height=50, units="mm")
# 


# Name <- gsub(" ", "_", OM@Name)
# df$Name <- OM@Name
# flname <- paste0("VmaxDF_", Name, ".rdata")
# saveRDS(df, file.path(Resultsdir, flname))
# 
# 
# for (Dep in DepVec) {
#   for (om in 1:length(OMlist)) {
#     runOM <- OMlist[[om]]
#     Hist <- runMSE(runOM, Hist=TRUE)
#     dep <- round(Hist$MSYs$SSBMSY_SSB0[1],3) * Dep 
#     runOM@D <- rep(dep, 2)
#     
#     name <- names(OMlist)[om]
#     message(OM@Name, " - ", name, " - Dep: ", Dep)
#     
#     MSE <- runMSE(runOM, MPs=testMPs, parallel = TRUE)
#     
#     omdat <- MSE@OM
#     omdat$sim <- 1:MSE@nsim
#     obsdat <-MSE@Obs
#     obsdat$sim <- 1:MSE@nsim
#     
#     output <- data.frame(B_BMSY=as.vector(unlist(MSE@B_BMSY)),
#                          F_FMSY=as.vector(unlist(MSE@F_FMSY)),
#                          Catch=as.vector(unlist(MSE@C)),
#                          TAC=as.vector(unlist(MSE@TAC)),
#                          MP=rep(MSE@MPs, each=MSE@nsim),
#                          sim=1:MSE@nsim, Years=rep(1:MSE@proyears, each=MSE@nsim*MSE@nMPs))
#     
#     outDF <- left_join(output, omdat, by='sim')
#     outDF <- left_join(outDF, obsdat, by='sim')
#     outDF$dep <- dep
#     outDF$Name <- MSE@Name
#     Name <- gsub(" ", "_", OM@Name)
#     flname <- paste0("Vmaxage_", Name, "_", name, "_", Dep, ".rdata")
#     saveRDS(outDF, file.path(Resultsdir, flname))
#   }
# }
# 
# 
