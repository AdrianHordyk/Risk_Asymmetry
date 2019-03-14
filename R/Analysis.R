
# devtools::install_github('tcarruth/MSEtool', ref = "new_dome")

source("R/control.r")
library(cowplot)

Y1 <- 4 
Y2 <- 30
lsize <- 0.5
lsize2 <- 0.3
cols1 <- c('#357DD9', "#C7FA2C")
cols2 <- c('#2B4970','#87A236')
text.size <- 12


# ---- Plot Maturity & Selectivity Curves ----
templist <- list()
for (x in 1:4) {
  OM <- MakeOM(x, nsim) # load OM 
  OM@K <- OM@K / 12 
  OM@t0 <- OM@t0 * 12 
  OM@maxage <- OM@maxage * 12 
  OM@M <- OM@M / 12
  OM@nsim <- 2
  Hist <- runMSE(OM, Hist=TRUE, silent = TRUE)
  Name <- OM@Name
  Length <- Hist@AtAge$Length[1,,1]
  Select <- Hist@AtAge$Select[1,,1]
  Maturity <- Hist@AtAge$Maturity[1,,1]
  Age <- 1:OM@maxage
  DF <- data.frame(Name=Name, Age=Age, Length=Length, Selectivity=Select, Maturity=Maturity)
  templist[[x]] <- DF %>% tidyr::gather("key", "value", 4:5)
}
DF <- do.call("rbind", templist)
DF$Name <- as.character(DF$Name)
DF <- DF %>% arrange((Name))
DF$Name <- as.factor(DF$Name)
DF$Age <- DF$Age/12
t1 <- DF %>% group_by(Name) %>% filter(key=="Maturity") %>%
  summarize(A50=Age[min(which(value>=0.5))],
            x=quantile(Age, 0.9),
            y=0.75)

t2 <- DF %>% group_by(Name) %>% filter(key=="Selectivity") %>%
  summarize(S50=Age[min(which(value>=0.5))])
tdf <- left_join(t1, t2) %>% mutate_at(c(2,5), round, 1)
tdf$key <- NA

DF %>% filter(Name=="Pacific ocean perch", key=='Maturity') 

p1 <- ggplot(DF, aes(x=Age, y=value,linetype=key)) + facet_wrap(~Name, scale="free") +
    geom_line(size=1) + labs(x="Age (year)", y="Proportion", linetype='') + 
  theme_classic() + theme(strip.background = element_blank()) +
  geom_text(data=tdf, aes(x=x, y=y, label=paste("A50 = ", A50)), size=3.5) +
  geom_text(data=tdf, aes(x=x, y=y*0.9, label=paste("S50 = ", S50)), size=3.5)
  
p1 

ggsave("Figures/Figure1.png", p1, units='mm', width=180, height=170, dpi=600)

# ---- Figure 1 ----
OM <- new("OM")
OM <- testOM
OM@maxage <- 40
OM@M <- c(0.2, 0.2)
OM@h <- c(0.7, 0.7)
OM@Linf <- c(100,100)
OM@K <- c(0.1, 0.1)
OM@L50 <- c(60, 60)
OM@L50_95 <- c(2,2)
OM@L5 <- c(55, 55)
OM@LFS <- c(58, 58)
OM@Vmaxlen <- c(1,1)
OM@isRel <- "FALSE"
OM@D <- c(0.2, 0.8)
OM <- tinyErr(OM)

OM@a <- 1E-5 
OM@b <- 3 

msys <- CalcMSYs(OM)
msys

YieldC <- CalcYieldCurve(OM) %>% mutate(Yield=Yield/msys$Yield)

Mbias <- exp(c(-0.25, 0.25))
msys_2 <- CalcMSYs(OM, M=OM@M * Mbias[1])
msys_3 <- CalcMSYs(OM, M=OM@M * Mbias[2])

msys
msys_2
msys_3


NegB <- EqCalcs(log(msys_2$F), OM,opt=2)  %>% t() %>% as.data.frame() %>%
  mutate(Yield=Yield/msys$Yield) %>% select(Yield, F, SB_SB0)

PosB <- EqCalcs(log(msys_3$F), OM,opt=2)  %>% t() %>% as.data.frame() %>%
  mutate(Yield=Yield/msys$Yield) %>% select(Yield, F, SB_SB0)

df <- data.frame(B=c(msys$SB_SB0, NegB$SB_SB0, PosB$SB_SB0),
                 F=c(msys$F, msys_2$F, msys_3$F),
                 Y=c(1, NegB$Yield, PosB$Yield),         
                 bias=c("Unbiased", "-25%", "+25%"))

lsize <- 1
size <- 3

Y_F <- ggplot(YieldC, aes(x=F, y=Yield)) + geom_line(size=lsize) +
  theme_classic() +
  labs(x="Fishing Mortality", y="Yield", shape="Bias log(M)")
Y_Fa <- Y_F +  geom_point(data=df, aes(x=F, y=Y, shape=fct_inorder(bias)), size=size) +
  geom_text(aes(x=0, y=1, label="a)"))

Y_D <- ggplot(YieldC, aes(x=SB_SB0, y=Yield)) + geom_line(size=lsize) + 
  theme_classic() + 
  labs(x="Spawning Depletion", y="Yield",shape="Bias log(M)")
Y_Da <- Y_D + geom_point(data=df, aes(x=B, y=Y, shape=fct_inorder(bias)), size=size) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  geom_text(aes(x=0, y=1, label="b)"))

join_plots(list(Y_Fa, Y_Da))

# ggsave("Figures/Presentation/Y_F1.png", PlotList$Y_F + theme, width=width, height=height)



# ---- Mbias Results ----

# --- Example Results for Fig 3 & 4 ----
Dep <- 0.5 
Sp <- "Pacific_hake"

makeMbiasDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("Mbias", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  Vars <- Data$Mbias %>% unique()
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(Mbias == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, Mbias)
    df$bias <- round(log(Vars[x]),3)
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(log(Vars[x])<0, "Negative", "Positive")
    if (log(Vars[x])==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF$BMSY <- Data$SSBMSY_SSB0 %>% unique()
  DF
}

DF <- makeMbiasDF('Pacific_hake')


tt <- DF %>% group_by(Mbias) %>% summarize(unique(Mbias))
DF %>% filter(Years > 10, Mbias==1) %>% summarize(median(D))

DF %>% filter(Years > 5, Mbias==tt$Mbias[1]) %>% summarize(median(D))
DF %>% filter(Years > 5, Mbias==tt$Mbias[9]) %>% summarize(median(D))

DF %>% filter(Years ==1, Mbias==tt$Mbias[1]) %>% summarize(median(D))
DF %>% filter(Years ==1, Mbias==1) %>% summarize(median(D))

DF %>% filter(Years ==30, Mbias>1) %>% group_by(Mbias) %>% summarize(median(stRelC))
DF %>% filter(Years ==30, Mbias<1) %>% group_by(Mbias) %>% summarize(median(stRelC))

DF %>% filter(Years ==30, Mbias>1) %>% group_by(Mbias) %>% summarize(median(relD))
DF %>% filter(Years ==30, Mbias<1) %>% group_by(Mbias) %>% summarize(median(relD))

dat <- DF %>% group_by(Years, bias) %>%
  summarise(Catch=median(RelC),
            Depletion=median(D)) %>%
  tidyr::gather("key", "value", 3:4)

dat$bias %>% as.factor %>% levels() 
linedf <- data.frame(key=c("Catch", 'Depletion'), yintercept=c(1, unique(round(DF$BMSY,2))))

cols <- gplots::rich.colors(9)

# linetype=as.factor(bias))
P1 <- ggplot(dat, aes(x=Years, y=value, color=as.factor(bias))) + 
  geom_hline(data=linedf, aes(yintercept=yintercept), linetype=4, col="darkgray", size=lsize2) +
  geom_line(size=lsize) +
  scale_color_manual(values=cols) +
  facet_wrap(~key, scales="free_y") +
  expand_limits(y=c(0,1)) +
  labs(y="Median", color="bias log(M)", linetype="bias log(M)",
       tag = "a)") + 
  guides(color=FALSE, linetype=FALSE) + theme_classic() +
  theme(strip.background = element_blank(),
        plot.tag = element_text(size=8))

dat <- DF %>% group_by(Years, bias) %>%
  summarise(Catch=median(stRelC),
            Depletion=median(relD)) %>%
  tidyr::gather("key", "value", 3:4)

P2 <- ggplot(dat, aes(x=Years, y=value, color=as.factor(bias))) + 
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray", size=lsize2) +
  geom_hline(yintercept = 2, linetype=4, col="darkgray", size=lsize2) +
  geom_line(size=lsize) +
  scale_color_manual(values=cols) +
  facet_wrap(~key, scales="free_y") +
  expand_limits(y=c(0,1)) +
  labs(y="Median", color="bias log(M)", linetype="bias log(M)",
       tag = "b)") + 
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.8, 'lines'),
        plot.tag = element_text(size=8))

leg <- cowplot::get_legend(P2)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust = -1)

Pout <- cowplot::plot_grid(Plots, leg, rel_widths = c(2, 0.2))

Pout

ggsave("Figures/Figure3.png", Pout, units='mm', width=180, height=170, dpi=600)


# ---- Figure 4 -----

dat <- DF %>% filter(Years %in% c(Y1, Y2)) %>% 
  select(sim, bias, Years, Catch=stRelC, Depletion=relD) %>% 
  tidyr::gather("key", "value", 4:5)

P3a <- ggplot(dat %>% filter(key=="Catch"),
              aes(x=as.factor(bias), y=value, fill=as.factor(Years))) +
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 4), expand=c(0,0)) +
  theme_classic() + theme(strip.background = element_blank()) + 
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90))+
  labs(fill="Year", x="Bias log(M)", y="Relative Change in Catch", tag="a)") +
  guides(fill=FALSE) +
  scale_fill_manual(values=cols1)

P3b <- ggplot(dat %>% filter(key!="Catch"),
              aes(x=as.factor(bias), y=value, fill=as.factor(Years))) + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 4), expand=c(0,0)) +
  theme_classic() + theme(strip.background = element_blank()) + 
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90))+
  labs(fill="Year", x="Bias log(M)", y="Relative Change in Depletion", tag="b)") +
  scale_fill_manual(values=cols1) 

P3 <- cowplot::plot_grid(P3a, P3b, ncol=2, rel_widths = c(0.5, 0.5))
title <- ggdraw() + 
  draw_label("Change in Catch and Spawning Depletion with bias in log(M)")
P3 <- plot_grid(title, P3, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

dat <- DF %>% group_by(abs_bias, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 4:5)

P4a <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
              aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 6), expand=c(0,0)) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90))+
  labs(fill="Year", x="Absolute bias log(M)", y="Catch Ratio", tag="c)") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE) 

P4b <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
              aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) +
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 2), expand=c(0,0)) +
  theme_classic() + theme(strip.background = element_blank()) + 
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90))+
  labs(fill="Year", x="Absolute bias log(M)", y="Depletion Ratio", tag="d)") +
  scale_fill_manual(values=cols2) 
 
P4 <- cowplot::plot_grid(P4a, P4b, ncol=2, rel_widths = c(0.5, 0.5))

title <- ggdraw() + 
  draw_label("Asymmetry in Catch and Spawning Depletion with over-estimation of log(M)")
P4 <- plot_grid(title, P4, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

Pout <- cowplot::plot_grid(P3, P4, nrow=2)

ggsave("Figures/Figure4_unedited.png", Pout, dpi=600, units='mm', width=180,
       height=170)  

# --- Compare all methods ----

funs <- list(makeMbiasDF, makehbiasDF, makeageMsDF, makebetasDF, makeCatchDF, makeSeletDF)
methods <- c("M", "CR", 'A50', 'Index', 'Catch', 'Vmax')
templist <- list()
for (x in seq_along(funs)) {
  Phake <- funs[[x]]('Pacific_hake')
  Arrow <- funs[[x]]('Arrowtooth_Flounder')
  Swarehou <- funs[[x]]('Silver_warehou')
  POP <- funs[[x]]('Pacific_ocean_perch')
  
  DF <- bind_rows(Phake, Arrow, Swarehou, POP)
  DF$Species <- gsub("_", " ", DF$Species)
  DF$Species <- tolower(DF$Species)
  DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))
  
  # Expected Loss with over-estimation
  dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
    summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
              Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
    filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)
  dat$Analysis <- methods[x]
  templist[[x]] <- dat
}
DF <- do.call('rbind', templist)
DF$Analysis <- factor(DF$Analysis, levels=methods, ordered = TRUE)

df <- DF %>% filter(abs_bias!=0) %>% group_by(Years, key, Analysis, Species, abs_bias) %>%
  summarize(med=median(value), low=quantile(value, 0.05), high=quantile(value, 0.95))

df$Analysis <- factor(df$Analysis, levels=methods, ordered = TRUE)

abs_bias <- unique(df$abs_bias)

linedf_a <- data.frame(key=c(rep("Catch",5), rep("Depletion", 4)),
                       yintercept=c(1, 0.5, 2, 0.25, 4, 1, 0.5, 2, 0.25),
                       linetype=c(3,4,4,5,5, 3,4,4,5),
                       col="darkgray")

linedf <- do.call("rbind", rep(list(linedf_a), length(methods)))
linedf$Analysis <- rep(methods, each=nrow(linedf_a))
linedf$Analysis <- factor(linedf$Analysis, levels=methods, ordered = TRUE)

Pout <- ggplot(df, aes(x=abs_bias, y=med, color=Species, linetype=as.factor(Years))) + 
  geom_line() + geom_point() + 
  geom_hline(data=linedf, aes(yintercept=yintercept), col="darkgray", linetype=linedf$linetype) +
  facet_grid(key~Analysis, scales = 'free_y') + 
  theme_classic() + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous("Absolute bias", labels = as.character(abs_bias), breaks = abs_bias) +
  expand_limits(y=0) + 
  coord_cartesian(clip = 'off') +
  labs(x="Absolute bias", y="Median Ratio", linetype='Year', color="Stock") +
  theme(strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"),
        legend.position = 'right',
        axis.text.x = element_text(angle=90))
Pout

ggsave("Figures/Figure5.png", Pout, dpi=600, width=180, height=120, units="mm")


# --- Mbias Results ----

Phake <- makeMbiasDF('Pacific_hake')
Arrow <- makeMbiasDF('Arrowtooth_Flounder')
Swarehou <- makeMbiasDF('Silver_warehou')
POP <- makeMbiasDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))


# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.99))
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.99))
ymax2$max <- max(ymax2$max, 1)
scaleFUN <- function(x) sprintf("%.2f", x)

# Relative change 
text.size <- 10
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
              aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray", size=1.1) +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text=element_text(size=text.size),
        strip.text.y=element_blank(),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Absolute bias log(M)", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray", size=1.1) +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax2$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y=element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Absolute bias log(M)", y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)


leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)


Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure6.png", Pout, dpi=600, width=180, height=160, units="mm")


tab <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))%>% 
  group_by(abs_bias, Species, Years) %>%
  summarize(medD = median(Depletion), D5 = quantile(Depletion, 0.05),
            D95 = quantile(Depletion, 0.95),
            medC = median(Catch), C5 = quantile(Catch, 0.05),
            C95 = quantile(Catch, 0.95))

tab %>% filter(Years==4, abs_bias>=0.375) %>% group_by(abs_bias, Species) %>% 
  summarize(min=min(medD), max=max(medD))
            
tab %>% filter(Years==4, Species=="Pacific ocean perch") %>% group_by(abs_bias, Species) %>% 
  summarize(min=min(medD), max=max(medD))
            

tab %>% filter(Years==30) %>% group_by(abs_bias, Species) %>% 
  summarize(min=min(medD), max=max(medD)) %>% 
  filter(Species=="Pacific ocean perch")


tab %>% filter(abs_bias==0.125,Years==4) %>% mutate_at(3:9, round, 2)
tab %>% filter(abs_bias==0.125, Years==30) %>% mutate_at(3:9, round, 2) 

tab %>% filter(Years==30) %>% mutate_at(3:9, round, 2) %>% 
  filter(abs_bias==.5) %>% select(medD) %>% 
  ungroup() %>% 
  summarize(m=min(medD), max(medD))


# ---- Hbias Results ----


makehbiasDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("Hbias", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  
  Vars <- Data$var %>% unique()
  
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(var == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, var)
    
    df$bias <- h2CR(Data$hs[1] + Vars[x])/h2CR(unique(Data$hs)) -1 
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(df$bias[1]<0, "Negative", "Positive")
    if (df$bias[1]==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF
}


Phake <- makehbiasDF('Pacific_hake')
Arrow <- makehbiasDF('Arrowtooth_Flounder')
Swarehou <- makehbiasDF('Silver_warehou')
POP <- makehbiasDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))

# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.99))
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.99))
ymax2$max <- max(ymax2$max, 1)


# Relative change 
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Bias CR", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax2$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Absolute bias CR", y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)


leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)

Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure7.png", Pout, dpi=600, width=180, height=160, units="mm")





tab <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))%>% 
  group_by(abs_bias, Species, Years) %>%
  summarize(medD = median(Depletion), D5 = quantile(Depletion, 0.05),
            D95 = quantile(Depletion, 0.95),
            medC = median(Catch), C5 = quantile(Catch, 0.05),
            C95 = quantile(Catch, 0.95))

tab %>% group_by(Species, Years) %>% 
  filter(Species=="Pacific hake", Years==30) %>% select(medD, medC)

tab %>% group_by(Species, Years) %>%
  summarize(min=min(medD), max=max(medD)) %>% 
  filter(Species=="Pacific ocean perch")


tab %>% filter(abs_bias==0.125,Years==4) %>% mutate_at(3:9, round, 1)

tab %>% filter(abs_bias==0.125, Years==30) %>% mutate_at(3:9, round, 1) 

tab %>% filter(Years==30) %>% mutate_at(3:9, round, 2) %>% 
  filter(abs_bias==.5) %>% select(medD) %>% 
  ungroup() %>% 
  summarize(m=min(medD), max(medD))



# --- Age of Maturity -----

makeageMsDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("AgeMbias", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  Vars <- Data$var %>% unique()
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(var == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, var)
    
    df$bias <- df$var - 1 
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(df$bias[1]<0, "Negative", "Positive")
    if (df$bias[1]==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF
}


Phake <- makeageMsDF('Pacific_hake')
Arrow <- makeageMsDF('Arrowtooth_Flounder')
Swarehou <- makeageMsDF('Silver_warehou')
POP <- makeageMsDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))

# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2$max <- max(ymax2$max, 1)


# Relative change 
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Bias A50", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax2$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="Absolute bias A50", y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)


leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)

Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure8.png", Pout, dpi=600, width=180, height=160, units="mm")



tab <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))%>% 
  group_by(abs_bias, Species, Years) %>%
  summarize(medD = median(Depletion), D5 = quantile(Depletion, 0.05),
            D95 = quantile(Depletion, 0.95),
            medC = median(Catch), C5 = quantile(Catch, 0.05),
            C95 = quantile(Catch, 0.95))

t1 <- tab %>% filter(Species =='Silver warehou', Years==4) %>% select(medC)
1 - t1$medC

tab %>% filter(abs_bias==0.125,Years==4) %>% mutate_at(3:9, round, 2)

tab %>% filter(abs_bias==0.125, Years==30) %>% mutate_at(3:9, round, 2) 

tab %>% group_by(Species, abs_bias) %>% filter(Years==4) %>%
  summarize(min=min(medC), max=max(medC)) %>% filter(Species == "Silver warehou") %>%
  group_by(abs_bias) %>% summarize(t=1-min)

tab %>% group_by(Species, Years) %>% filter(abs_bias==0.5) %>%
  summarize(min=min(medD), max=max(medD)) %>% 
  filter(Species=="Pacific ocean perch")


temp <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))

temp <- temp %>% filter(abs_bias==0.125, Years==30, Species=="Silver warehou")
hist(temp$Catch)
quantile(temp$Catch)


tab %>% filter(Years==30) %>% mutate_at(3:9, round, 2) %>% 
  filter(abs_bias==.5) %>% select(medD) %>% 
  ungroup() %>% 
  summarize(m=min(medD), max(medD))

# ---- Index of Abundance ---- 

makebetasDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("betabias", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  Vars <- Data$var %>% unique()
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(var == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, var)
    
    df$bias <- log(df$var)
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(df$bias[1]<0, "Negative", "Positive")
    if (df$bias[1]==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF
}


Phake <- makebetasDF('Pacific_hake')
Arrow <- makebetasDF('Arrowtooth_Flounder')
Swarehou <- makebetasDF('Silver_warehou')
POP <- makebetasDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))

# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2$max <- max(ymax2$max, 1)


# Relative change 
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 2), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="\beta", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 4), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x=expression('Absolute log ' * beta), y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)

leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)

Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure9.png", Pout, dpi=600, width=180, height=160, units="mm")


tab <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))%>% 
  group_by(abs_bias, Species, Years) %>%
  summarize(medD = median(Depletion), D5 = quantile(Depletion, 0.05),
            D95 = quantile(Depletion, 0.95),
            medC = median(Catch), C5 = quantile(Catch, 0.05),
            C95 = quantile(Catch, 0.95))


tab %>% filter(Years==4, Species=="Pacific ocean perch") %>%
  select(medC)


# --- Reporting of Catch ----

makeCatchDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("Cbias", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  Vars <- Data$var %>% unique()
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(var == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, var)
    
    df$bias <- 1-(df$var)
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(df$bias[1]<0, "Negative", "Positive")
    if (df$bias[1]==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF
}


Phake <- makeCatchDF('Pacific_hake')
Arrow <- makeCatchDF('Arrowtooth_Flounder')
Swarehou <- makeCatchDF('Silver_warehou')
POP <- makeCatchDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))

# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2$max <- max(ymax2$max, 1)


# Relative change 
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="\beta", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax2$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x=paste('Absolute catch bias'), y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)

leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)

Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure10.png", Pout, dpi=600, width=180, height=160, units="mm")


tab <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  filter(Years%in%c(4,30), abs_bias!=0) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]]))%>% 
  group_by(abs_bias, Species, Years) %>%
  summarize(medD = median(Depletion), D5 = quantile(Depletion, 0.05),
            D95 = quantile(Depletion, 0.95),
            medC = median(Catch), C5 = quantile(Catch, 0.05),
            C95 = quantile(Catch, 0.95))


tab %>% filter(Years==4, Species=="Pacific ocean perch")



# ---- Selectivity ----

makeSeletDF <- function(Sp, Dep=0.5) {
  Data <- readRDS(file.path("Results", paste0(paste("select", Sp, Dep, sep="_"), '.rdata')))
  Data <- Data %>% mutate(D=((B_BMSY * SSBMSY) / SSB0), RelC=Catch/RefY)
  Data$MP <- as.character(Data$MP)
  Vars <- Data$MP %>% unique() 
  List <- list(); df <- NA
  for (x in seq_along(Vars)) {
    rm(df)
    df <- Data %>% dplyr::filter(MP == Vars[x]) %>% select(sim, D, B_BMSY, RelC, Years, MP)
    val <- as.numeric(unlist(regmatches(Vars[x],
                                        gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                 Vars[x]))))
    df$bias <- val / 0.65 - 1
    df$abs_bias <- round(abs(df$bias), 3)
    df$Bias <- ifelse(df$bias[1]<0, "Negative", "Positive")
    if (df$bias[1]==0) {
      df2 <- df
      df2$Bias <- "Negative"
      df$Bias <- "Positive"
      df <- bind_rows(df2, df)
    }
    List[[x]] <- df
  }
  DF <- do.call("rbind", List)
  DF$Species <- Sp
  DF$relD <- DF$D/DF$D[DF$abs_bias==0]
  DF$stRelC <- DF$RelC/DF$RelC[DF$abs_bias==0]
  DF
}


Phake <- makeSeletDF('Pacific_hake')
Arrow <- makeSeletDF('Arrowtooth_Flounder')
Swarehou <- makeSeletDF('Silver_warehou')
POP <- makeSeletDF('Pacific_ocean_perch')

DF <- bind_rows(Phake, Arrow, Swarehou, POP)
DF$Species <- gsub("_", " ", DF$Species)
DF$Species <- tolower(DF$Species)
DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))

# Expected Loss with over-estimation
dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
  summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
            Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
  filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)

ymax <- dat %>% filter(key=="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax$max <- 1.5
ymax2 <- dat %>% filter(key!="Catch") %>% ungroup() %>% summarize(max=quantile(value, 0.95))
ymax2$max <- max(ymax2$max, 1)


# Relative change 
P1 <- ggplot(dat %>% filter(key=="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        legend.position = 'top',
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x="", y="Catch Ratio") +
  scale_fill_manual(values=cols2) 

P2 <- ggplot(dat %>% filter(key!="Catch", abs_bias>0), 
             aes(x=as.factor(abs_bias), y=value, fill=as.factor(Years))) + 
  facet_grid(key~Species, scales='free_y') + 
  geom_hline(yintercept = 1, linetype=3, col="darkgray") +
  geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
  geom_hline(yintercept = 2, linetype=4, col="darkgray") +
  geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
  geom_hline(yintercept = 4, linetype=5, col="darkgray") +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, ymax2$max), expand=c(0,0),labels=scaleFUN) +
  theme_classic() + theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size=text.size),
        axis.text = element_text(size=text.size-1),
        axis.text.x = element_text(angle=90),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.text=element_text(size=text.size),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(fill="Year", x=paste('Absolute bias Vmax'), y="Depletion Ratio") +
  scale_fill_manual(values=cols2) + guides(fill=FALSE)

leg <- cowplot::get_legend(P1)
Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                            P2 + theme(legend.position="none"),
                            align = 'vh',
                            ncol=1, nrow=2,
                            hjust =-1)

Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
ggsave("Figures/Figure11.png", Pout, dpi=600, width=180, height=160, units="mm")



# Supp Figures ####
Y1 <- 4 ; Y2 <- 30 

funs <- list(makeMbiasDF, makehbiasDF, makeageMsDF, makebetasDF, makeCatchDF, makeSeletDF)
methods <- c("M", "CR", 'A50', 'Index', 'Catch', 'Vmax')
xlab <- c("Bias log(M)",
          "Bias CR",
          'Bias A50',
          'Bias log(beta)',
          'Bias Catch reporting',
          'Bias Vmax')
text.size <- 8 
for (x in seq_along(funs)) {
  Phake <- funs[[x]]('Pacific_hake')
  Arrow <- funs[[x]]('Arrowtooth_Flounder')
  Swarehou <- funs[[x]]('Silver_warehou')
  POP <- funs[[x]]('Pacific_ocean_perch')
  
  DF <- bind_rows(Phake, Arrow, Swarehou, POP)
  DF$Species <- gsub("_", " ", DF$Species)
  DF$Species <- tolower(DF$Species)
  DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))
  
  dat <- DF %>% filter(Years %in% c(Y1, Y2)) %>% 
    select(sim, bias, Years, Catch=stRelC, Depletion=relD, Species) %>% 
    tidyr::gather("key", "value", 4:5)
  dat$bias <- round(dat$bias, 2)
  Pout <- ggplot(dat,
                aes(x=as.factor(bias), y=value, fill=as.factor(Years))) +
    facet_grid(key~Species, scales = 'free') +
    geom_hline(yintercept = 1, linetype=3, col="darkgray") +
    geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
    geom_hline(yintercept = 2, linetype=4, col="darkgray") +
    geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
    geom_hline(yintercept = 4, linetype=5, col="darkgray") +
    geom_boxplot(outlier.shape = NA) + 
    scale_y_continuous(limits = c(0, 4), expand=c(0,0)) +
    theme_classic() + theme(strip.background = element_blank()) + 
    theme(axis.title = element_text(size=text.size),
          axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90))+
    labs(fill="Year", x=xlab[x], y="Relative Change") +
    scale_fill_manual(values=cols1)
  
  filename <- paste0("Figures/SFigure", x, ".png")
  ggsave(Pout, filename = filename, dpi=600, width=180, height=160, units="mm")
}



# ----------------- END OF ANALYSIS -------------------------------




