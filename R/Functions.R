
#### Generate Operating Models ####
MakeOM <- function(omnum, nsim=100, interval=4, nyears=50, proyears=30, obs_mod=Precise_Unbiased) {
  dat <- readxl::read_xlsx("OMs/OMs.xlsx")[omnum,]
  
  OM <- testOM # copy & modify 
  OM@nsim <- nsim
  OM@Species <- dat$Species
  OM <- tinyErr(OM, silent=TRUE)
  OM <- Replace(OM, obs_mod)
  OM@Name <- dat$Name
  
  OM@AC <- c(0,0)
  OM@isRel <- "FALSE"
  nms <- names(dat)
  nms <- nms[!nms %in% c('Name', "Species")]
  for (nm in nms) {
    slot(OM, nm) <- rep(dat[[nm]],2)
  }
  
  OM@maxF <- 2
  OM@CurrentYr <- 2018
  
  OM@maxage <- ceiling(-log(0.01)/mean(OM@M))
  OM@nyears <- nyears
  OM@proyears <- proyears
  
  OM@EffYears <- c(1, ceiling(0.25*nyears), ceiling(0.75*nyears), nyears)
  OM@EffLower <- c(0, 0.75, 1, 1)
  OM@EffUpper <- c(0, 0.75, 1, 1)
  
  A50 <- (log(1 - (dat$L50/dat$Linf))/-dat$K) + dat$t0
  OM@SizeLimSD <- c(A50,A50) # used internally
  
  OM
  
}


matA <- function(A50, A95, ages) 1/(1+exp((-log(19)*(ages-A50)/(A95-A50))))

selA <- function(AS5, AFS, Vmaxage, maxage, nsim, OM) {
  
  srs <- (maxage - AFS) / ((-log(Vmaxage,2))^0.5) # 
  sls <- (AFS - AS5) /((-log(0.05,2))^0.5)
  srs <- rep(srs, nsim)
  sls <- rep(sls, nsim)
  AFS <- rep(AFS, nsim)
  
  vuln <- t(sapply(1:nsim, DLMtool:::getsel, lens=matrix(1:maxage, nrow=nsim, ncol=maxage, byrow=TRUE), 
                   lfs=AFS, sls=sls, srs=srs))
  V <- array(vuln, dim=c(nsim, maxage, OM@nyears + OM@proyears))
  V
}

VB <- function(Linf, K, t0, ages) Linf * (1-exp(-K*(ages-t0)))
  
Calct0 <- function(Linf, K, L_Amin, Amin) {
  opt <- optimise(function(t0) (VB(Linf, K, t0, Amin) - L_Amin)^2, interval=c(-5,5))
  opt$minimum
}

# PopOM <- function(Name="name", Species, M, Linf, K, t0=NULL, L_Amin, Amin, h, Perr, A50=NULL, A95=NULL,
#                   L50=NULL, L95=NULL, slope=NULL, AS5=NULL, AFS=NULL, L5=NULL, LFS=NULL, Vmaxlen, 
#                   LR5=NULL, LFR=NULL, Rmaxlen=NULL, nsim, proyears=50) {
#   OM <- new("OM")
#   OM@Name <- "blank"
#   OM@Species <- Species
#   
#   OM@Fdisc <- c(1,1)
#   OM@maxF <- 0.8
#   OM@reps <- 1
#   OM@R0 <- 1E4
#   OM@SRrel <- 1 
#   OM@CurrentYr <- 2018
#   OM@Spat_targ <- c(1,1)
#   OM <- Replace(OM, Perfect_Info, silent = TRUE)
#   OM <- Replace(OM, Perfect_Imp, silent=TRUE)
#   OM <- tinyErr(OM, silent=TRUE)
# 
#   OM@nsim <- nsim 
#   OM@interval <- 4
#   OM@pstar <- 0.5
#   OM@L50 <- OM@L50_95 <- OM@L5 <- OM@LFS <- OM@Vmaxlen <- c(0,0)
#   OM@D <- c(0.5,0.5) # temp value
#   OM@isRel <- 'FALSE'
#   OM@AC <- c(0,0)
#   OM@Frac_area_1 <- OM@Size_area_1 <- OM@Prob_staying <- c(0.5,0.5)
#   OM@a <- 1E-5 
#   OM@b <- 3
#   OM@Msd <- c(0,0)
#   
#   OM@Name <- Name
#   OM@M <- rep(M, 2)
#   OM@Linf <- rep(Linf, 2)
#   OM@K <- rep(K, 2)
#   if (is.null(t0)) {
#     t0 <- Calct0(Linf, K, L_Amin, Amin)
#     print(Name)
#     print(t0)
#   }
#   
#   OM@t0 <- rep(t0, 2)
#   OM@h <- rep(h, 2)
#   OM@Perr <- rep(Perr, 2)
#   OM@maxage <- maxage <- ceiling(-log(0.01)/mean(OM@M*0.75))
#   OM@nyears <- nyears <- ceiling(mean(OM@maxage) * 1.5)
#  
#   OM@proyears <- proyears # ceiling(mean(OM@maxage) * 2)
#   
#   OM@EffYears <- c(1, ceiling(0.25*nyears), ceiling(0.75*nyears), nyears)
#   OM@EffLower <- c(0, 0.75, 1, 1)
#   OM@EffUpper <- c(0, 0.75, 1, 1)
#   
#   ages <- 1:maxage
#   
#   lenage <- Linf * (1-exp(-K*(ages-t0)))
#   
#   if (is.null(A50)) {
#     A50 <- (log(1 - (L50/Linf))/-K) + t0
#     if (is.null(L95)) {
#       A95 <- log(19)/-slope + A50
#     } else {
#       A95 <- (log(1 - (L95/Linf))/-K) + t0  
#     } 
#     print(Name)
#     print(data.frame(A50=A50, A95=A95))
#   }
#   if (is.null(AS5)) {
#     AS5 <- (log(1 - (L5/Linf))/-K) + t0
#     AFS <- (log(1 - (LFS/Linf))/-K) + t0
#     print(Name)
#     print(data.frame(AS5=AS5, AFS=AFS))
#   }
#   
#   if (!is.null(LR5)) {
#     OM@LR5 <- rep(LR5, 2)
#     OM@LFR <- rep(LFR, 2)
#     OM@Rmaxlen <- rep(Rmaxlen, 2)
#   }
#   
#   OM@cpars$Mat_age <- aperm(array(matA(A50, A95, ages), dim=c(OM@maxage, OM@nsim, OM@nyears+OM@proyears)), c(2,1,3))
#   OM@cpars$V <- selA(AS5, AFS, Vmaxlen, maxage, nsim, OM)
#   OM@SizeLimSD <- c(A50,A50)
#   OM
#   
# }



EqCalcs <- function(logapicF, OM, Hist=NULL, opt=1, sim=1) {
  # Box 3.1 Walters & Martell 2004
  if (is.null(Hist)) {
    OM@nsim <- 2 
    Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  }
  maxage <- OM@maxage
  SRrel <- OM@SRrel
  M_at_Age <- Hist$SampPars$M_ageArray[sim,,OM@nyears]
  Len_at_Age <- Hist$SampPars$Len_age[sim,,OM@nyears]
  Wght_at_Age <- Hist$SampPars$Wt_age[sim,,OM@nyears]
  V_at_Age <- Hist$SampPars$V[sim,,OM@nyears]
  Mat_at_Age <- Hist$SampPars$Mat_age[sim,,OM@nyears]
  h <- Hist$SampPars$hs[sim]
  M <- max(M_at_Age)
  R0 <- OM@R0[1]
  
  apicF <- exp(logapicF)
  lx <- l0 <- rep(1, maxage)
  for (a in 2:maxage) {
    l0[a] <- l0[a-1] * exp(-M_at_Age[a-1])
    lx[a] <- lx[a-1] * exp(-(M_at_Age[a-1] + apicF*V_at_Age[a-1]))
  }
  Egg0 <- sum(l0 * Wght_at_Age * Mat_at_Age) # unfished egg production (assuming fecundity proportional to weight)
  EggF <- sum(lx * Wght_at_Age * Mat_at_Age) # fished egg production (assuming fecundity proportional to weight)
  
  vB0 <- sum(l0 * Wght_at_Age * V_at_Age)
  vBF <- sum(lx * Wght_at_Age * V_at_Age)
  
  SB0 <- sum(l0 * Wght_at_Age * Mat_at_Age) # same as eggs atm
  SBF <- sum(lx * Wght_at_Age * Mat_at_Age)
  
  B0 <- sum(l0 * Wght_at_Age) # same as eggs atm
  BF <- sum(lx * Wght_at_Age)
  
  Fage <- apicF * V_at_Age 
  Zage <- Fage + M_at_Age 
  
  YPR <-sum(Fage/Zage * (lx * Wght_at_Age ) * (1-exp(-Zage)))
  # Cata <- Fs/Zs * Wght_at_Age * N[,yr] * (1-exp(-Zs))
  
  h[h>0.999] <- 0.999
  recK <- (4*h)/(1-h) # Goodyear compensation ratio
  reca <- recK/Egg0
  if (SRrel ==1) recb <- (reca * Egg0 - 1)/(R0*Egg0) # BH SRR
  if (SRrel ==2) recb <- log(reca*Egg0)/(R0*Egg0) # Ricker SRR
  
  RelRec <- (reca * EggF-1)/(recb*EggF)
  RelRec[RelRec<0] <- 0
  
  Yield <- YPR * RelRec
  
  if (opt == 1)  return(-Yield)
  if (opt == 3)  return(Yield)
  if (opt == 2) {
    out <- c(Yield=Yield,
             F= -log(1 - (Yield/(vBF*RelRec*exp(-0.5*mean(M_at_Age))))),
             SB = SBF * RelRec,
             SB_SB0 = (SBF* RelRec)/(SB0*R0),
             B_B0 = BF/B0,
             B = BF * RelRec + Yield,
             VB = vBF * RelRec + Yield,
             VB_VB0 = vBF/vB0,
             RelRec=RelRec,
             SB0 = SB0 * R0)
    
    return(out)
  }
  
  
  return(NA)
}


PopModel <- function(OM, nyrs, apicF, type="all", sim=1) {
  OM@nsim <- 2 
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  maxage <- OM@maxage
  R0 <- OM@R0[1]
  hs <- OM@h[sim]
  SRrel <- OM@SRrel
  M_at_Age <- Hist$SampPars$M_ageArray[sim,,OM@nyears]
  Len_at_Age <- Hist$SampPars$Len_age[sim,,OM@nyears]
  Wght_at_Age <- Hist$SampPars$Wt_age[sim,,OM@nyears]
  V_at_Age <- Hist$SampPars$V[sim,,OM@nyears]
  Mat_at_Age <- Hist$SampPars$Mat_age[sim,,OM@nyears]
  
  if (length(apicF)==1) apicF <- rep(apicF, nyrs)
  N <- matrix(0, nrow=maxage, ncol=nyrs)
  SB <- Catch <- rep(0, nyrs)
  
  N[1,1] <- R0 
  N[2:maxage,1] <- N[1,1] * exp(-cumsum(M_at_Age[1:(maxage-1)]))
  SB[1] <- sum(N[,1] * Wght_at_Age * Mat_at_Age)
  SB0 <- SB[1]
  Egg0 <- SB[1]/R0
  hs[hs>0.999] <- 0.999
  
  recK <- (4*hs)/(1-hs) # Goodyear compensation ratio
  reca <- recK/Egg0
  if (SRrel ==1) recb <- (reca * Egg0 - 1)/(R0*Egg0) # BH SRR
  if (SRrel ==2) recb <- log(reca*Egg0)/(R0*Egg0) # Ricker SRR
  
  # Equilbrium 
  lx <- l0 <- rep(1, maxage)
  for (a in 2:maxage) {
    l0[a] <- l0[a-1] * exp(-M_at_Age[a-1])
    lx[a] <- lx[a-1] * exp(-(M_at_Age[a-1] + apicF[length(apicF)]*V_at_Age[a-1]))
  }
  
  EggF <- sum(lx * Wght_at_Age * Mat_at_Age) # fished egg production (assuming fecundity proportional to weight)
  RelRec <- (reca * EggF-1)/(recb*EggF)
  RelRec[RelRec<0] <- 0
  F <- rep(NA, nyrs)
  # Dynamic 
  for (yr in 2:nyrs) {
    N[1,yr] <- SRR(SB[yr-1], reca, recb, SRrel) 
    ind <- as.matrix(expand.grid(2:maxage, yr, 1:maxage, yr-1))
    for (a in 2:maxage) {
      N[a,yr] <- N[a-1,yr-1] * exp(-(M_at_Age[a-1]+V_at_Age[a-1]*apicF[yr]))
    }
    SB[yr] <- sum(N[,yr] * Wght_at_Age * Mat_at_Age)
    Fs <- V_at_Age*apicF[yr]
    Zs <- Fs + M_at_Age
    Cata <- Fs/Zs * Wght_at_Age * N[,yr] * (1-exp(-Zs))
    Catch[yr] <- sum(Cata)
    F[yr] <- -log(1-( Catch[yr]/(sum(N[,yr] * Wght_at_Age*V_at_Age *exp(-0.5*mean(M_at_Age))))))
  }
  
  
  if (type=='eq') return(RelRec)
  if (type=='dyn') return(N[1,nyrs])
  if (type =='all') return(data.frame(Eq=RelRec, Dyn=N[1,nyrs], 
                                      SB=SB[yr]/SB[1],
                                      Catch=Catch[yr]))
  if (type=='annual') 
    return(data.frame(Year=1:nyrs, Rec=N[1,], SB=SB, Catch=Catch,
                      F=F))
  
}

h2CR <- function(h)(4*h)/(1-h)
CR2h <- function(CR) CR/(CR+4)


CalcYieldCurve <- function(OM, sim=1) {
  OM@nsim <- 2 
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  
  M_at_Age <- Hist$SampPars$M_ageArray[sim,,OM@nyears]
  Len_at_Age <- Hist$SampPars$Len_age[sim,,OM@nyears]
  Wght_at_Age <- Hist$SampPars$Wt_age[sim,,OM@nyears]
  V_at_Age <- Hist$SampPars$V[sim,,OM@nyears]
  Mat_at_Age <- Hist$SampPars$Mat_age[sim,,OM@nyears]
  h <- Hist$SampPars$hs[sim]
  M <- max(M_at_Age)
  R0 <- OM@R0[1]
  # Calculate max F 
  MaxF <- 5 * M 
  
  Fvec <- seq(0, MaxF, length.out=100)
  Run <- sapply(seq_along(Fvec), function(x) 
    EqCalcs(logapicF=log(Fvec[x]), OM=OM, Hist=Hist, opt=2)) %>% t() %>% as.data.frame()
  

  
  Run <- Run %>% dplyr::filter(is.finite(Run$F)==TRUE)
  
  # Standardise Yield
  # Run <- Run %>% mutate(Yield=Yield/max(Yield))
  Run
}


CalcMSYs <- function(OM, sim=1, M=NULL, h=NULL) {
  OM@nsim <- 2
  
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  
  if (!is.null(h)) OM@h <- h
  if (!is.null(M)) OM@M <- M
  
  Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
  Hist$MSYs

   FMSY <-  Hist$MSYs$FMSY[1] # exp(optimise(EqCalcs, interval=log(c(0.01, 10*mean(OM@M))), 
  #                      OM=OM, Hist=Hist, opt=1)$minimum)
  
  EqCalcs(log(FMSY), OM=OM, Hist=Hist, opt=2) %>% 
    t() %>% as.data.frame()
  
}

SRR <- function(SB, reca, recb, SRrel) {
  if (SRrel == 1) return(reca * SB / (1+recb*SB)) # BHH
  if (SRrel == 2) return(reca * SB * exp(-recb*SB)) # Ricker
  return(NULL)
  
}


# https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(plots, layout_matrix, position = c("bottom", "right")) {
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  # gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight), layout_matrix=layout_matrix),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}


addLabels <- function(p, labelR=NULL, labelT=NULL, font.size=14, plot=FALSE) {
  
  z <- ggplotGrob(p)
  
  if (!is.null(labelR)) {
    # Get the positions of the strips in the gtable: t = top, l = left, ..
    posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
    # Add a new column to the right of current right strips, 
    # and a new row on top of current top strips
    width <- z$widths[max(posR$r)]    # width of current right strips
    z <- gtable_add_cols(z, width, max(posR$r))
    # Construct the new strip grobs
    stripR <- gTree(name = "Strip_right", children = gList(
      rectGrob(gp = gpar(col = NA, fill = "white")),
      textGrob(labelR, rot = -90, gp = gpar(fontsize = font.size, col = "grey10"))))
    # Position the grobs in the gtable
    z <- gtable_add_grob(z, stripR, t = min(posR$t), l = max(posR$r)+1, b = max(posR$b), name = "strip-right")
    
    # Add small gaps between strips
    z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  }
  if (!is.null(labelT)) {
    posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
    height <- z$heights[min(posT$t)]  # height of current top strips
    z <- gtable_add_rows(z, height, min(posT$t)-1)
    
    stripT <- gTree(name = "Strip_top", children = gList(
      rectGrob(gp = gpar(col = NA, fill = "white")),
      textGrob(labelT, gp = gpar(fontsize = font.size, col = "grey10"))))
    z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
    z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  }
  
  if (plot) grid.draw(z)
  return(invisible(z))
  
}

MakePlots <- function(OM, msys, msys_2, msys_3, size=2, lsize=1.25) {
  
  NegB <- EqCalcs(log(msys_2$F), OM,opt=2)  %>% t() %>% as.data.frame() %>%
    mutate(Yield=Yield/msys$Yield) %>% select(Yield, F, SB_SB0)
  
  PosB <- EqCalcs(log(msys_3$F), OM,opt=2)  %>% t() %>% as.data.frame() %>%
    mutate(Yield=Yield/msys$Yield) %>% select(Yield, F, SB_SB0)
  
  df <- data.frame(B=c(msys$SB_SB0, NegB$SB_SB0, PosB$SB_SB0),
                   F=c(msys$F, msys_2$F, msys_3$F),
                   Y=c(1, NegB$Yield, PosB$Yield),         
                   bias=c("Unbiased", "Negative", "Positive"))
  
  YieldC <- CalcYieldCurve(OM) %>% mutate(Yield=Yield/msys$Yield)
  
  Y_D <- ggplot(YieldC, aes(x=SB_SB0, y=Yield)) + geom_line(size=lsize) + 
    theme_classic() + 
    labs(x="Spawning Depletion", y="Yield",shape="Bias")
  Y_Da <- Y_D + geom_point(data=df, aes(x=B, y=Y, shape=fct_inorder(bias)), size=size) 
  
  Y_F <- ggplot(YieldC, aes(x=F, y=Yield)) + geom_line(size=lsize) +
    theme_classic() +
    labs(x="Fishing Mortality", y="Yield", shape="Bias")
  Y_Fa <- Y_F +  geom_point(data=df, aes(x=F, y=Y, shape=fct_inorder(bias)), size=size) 
  
  F_D <- ggplot(YieldC, aes(x=F, y=SB_SB0))+ geom_line(size=lsize) + 
    theme_classic() + 
    labs(x="Fishing Mortality", y="Spawning Depletion", shape="Bias")
  F_Da <- F_D + geom_point(data=df, aes(x=F, y=B, shape=fct_inorder(bias)), size=size) 
  
  # Projection plots 
  nyrs <- 100
  st <- 1.5
  Ftrend1 <- c(rep(st*msys$F, 0.5*nyrs), rep(msys$F,  0.5*nyrs))
  Ftrend2 <- c(rep(st*msys$F,  0.5*nyrs), rep(msys_2$F,  0.5*nyrs))
  Ftrend3 <- c(rep(st*msys$F,  0.5*nyrs), rep(msys_3$F,  0.5*nyrs))
  P1 <- PopModel(OM, nyrs, Ftrend1, type="annual")
  P1$bias <- "Unbiased"
  P2 <- PopModel(OM, nyrs, Ftrend2, type="annual")
  P2$bias <- "Negative"
  P3 <- PopModel(OM, nyrs, Ftrend3, type="annual")
  P3$bias <- "Positive"
  
  PDF <- bind_rows(P1, P2, P3)
  
  PDF <- PDF %>% mutate(Biomass=SB/max(SB), RelC=Catch/msys$Yield)
  
  df2 <- PDF %>% filter(Year==max(Year)) %>% select(RelC, Biomass, bias, Year)
  
  P_C <- ggplot(filter(PDF, Year>50), aes(x=Year, y=RelC, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
    expand_limits(y=0,1) + theme_classic() + 
    labs(x="Year", y="Yield", shape="Bias", linetype='Bias') +
    scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10))
  P_Ca <- P_C + geom_point(data=df2, aes(x=Year, y=RelC, shape=fct_inorder(bias)), size=size) 
  
  P_B <- ggplot(filter(PDF, Year>50), aes(x=Year, y=SB/msys$SB0, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
    expand_limits(y=c(0,1)) + theme_classic() +
    labs(x="Year", y="Spawning Depletion", shape="Bias", linetype='Bias') +
    scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10))
  P_Ba <- P_B + geom_point(data=df2, aes(x=Year, y=Biomass, shape=fct_inorder(bias)), size=size) 
  
  return(list(Y_F=Y_F, Y_Fa=Y_Fa,
              Y_D=Y_D, Y_Da=Y_Da,
              F_D=F_D, F_Da=F_Da,
              P_C=P_C, P_Ca=P_Ca,
              P_B=P_B, P_Ba=P_Ba,
              df=df, df2=df2,
              PDF=PDF))
}









