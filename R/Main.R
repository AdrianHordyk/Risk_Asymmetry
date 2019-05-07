
source("R/control.r")

Y1 <- 4 
Y2 <- 30
lsize <- 0.5
lsize2 <- 0.3
cols1 <- c('#357DD9', "#C7FA2C")
cols2 <- c('#2B4970','#87A236')
text.size <- 12

# ---- Figure 3 -----
DF <- makeMbiasDF('Pacific_hake')

dat <- DF %>% group_by(Years, bias) %>%
  summarise(Catch=median(RelC),
            Depletion=median(D)) %>%
  tidyr::gather("key", "value", 3:4)

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

ggsave("Figures/Figure3.png", Pout, units='mm', width=180, height=170, dpi=600)

dat %>% filter(key=="Catch", bias>0, Years==30)
dat %>% filter(key=="Catch", bias<0, Years==30)

dat %>% filter(key!="Catch", bias>0, Years==30)
dat %>% filter(key!="Catch", bias<0, Years==30)

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
  draw_label("Change in Catch and Depletion with bias in log(M)")
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
  draw_label("Asymmetry in Catch and Depletion with positive bias in log(M)")
P4 <- plot_grid(title, P4, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

Pout <- cowplot::plot_grid(P3, P4, nrow=2)

ggsave("Figures/Figure4.png", Pout, dpi=600, units='mm', width=180, height=170) 


# ---- Figures 5 & 6 - Medians for all Analyses ----

minVal <- 0.05 # values are truncated to this minimum value 

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

linedf <- do.call("rbind", rep(list(linedf_a), length(methods)*2))
linedf$Analysis <- rep(methods, each=nrow(linedf_a)*2)
linedf$Analysis <- factor(linedf$Analysis, levels=methods, ordered = TRUE)
linedf$Years <- rep(unique(df$Years), each=nrow(linedf_a))
linedf$Years <- factor(linedf$Years, ordered = TRUE)

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

df2 <- df
df2$med[df2$med < minVal] <- minVal

# --- Catch ---
df3 <- df2 %>% filter(key=='Catch')

df3$Years <- factor(df3$Years, ordered = TRUE)
Pout <- ggplot(df3, aes(x=abs_bias, y=med, color=Species)) + 
  geom_line() + geom_point() +
  facet_grid(Years~Analysis, scales = 'free_y') + 
  geom_hline(data=linedf, aes(yintercept=yintercept), col="darkgray", linetype=linedf$linetype) +
  theme_classic() + 
  scale_y_continuous(trans=log_trans(), breaks = base_breaks()) +
  # scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous("Absolute bias", labels = as.character(abs_bias), breaks = abs_bias) +
  expand_limits(y=minVal) + 
  coord_cartesian(clip = 'off') +
  labs(x="Absolute bias", y="Median Catch Ratio", linetype='Year', color="Stock") +
  theme(strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"),
        legend.position = 'right',
        axis.text.x = element_text(angle=90))

Pout <- addLabels(Pout, labelR="Year", labelT="")
ggsave("Figures/Figure5.png", Pout, dpi=600, width=180, height=120, units="mm")

# --- Depletion ---
df3 <- df2 %>% filter(key!='Catch')
df3$Years <- factor(df3$Years, ordered = TRUE)
Pout <- ggplot(df3, aes(x=abs_bias, y=med, color=Species)) + 
  geom_line() + geom_point() +
  facet_grid(Years~Analysis, scales = 'free_y') + 
  geom_hline(data=linedf, aes(yintercept=yintercept), col="darkgray", linetype=linedf$linetype) +
  theme_classic() + 
  scale_y_continuous(trans=log_trans(), breaks = base_breaks()) +
  # scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous("Absolute bias", labels = as.character(abs_bias), breaks = abs_bias) +
  expand_limits(y=minVal) + 
  coord_cartesian(clip = 'off') +
  labs(x="Absolute bias", y="Median Depletion Ratio", linetype='Year', color="Stock") +
  theme(strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"),
        legend.position = 'right',
        axis.text.x = element_text(angle=90))

Pout <- addLabels(Pout, labelR="Year", labelT="")
ggsave("Figures/Figure6.png", Pout, dpi=600, width=180, height=120, units="mm")



# ---- Supp Figures ----
funs <- list(makeMbiasDF, makehbiasDF, makeageMsDF, makebetasDF, makeCatchDF, makeSeletDF)
methods <- c("M", "CR", 'A50', 'Index', 'Catch', 'Vmax')
xlab <- c("Bias log(M)",
          "Bias CR",
          'Bias log(A50)',
          'Bias log(beta)',
          'Bias Catch reporting',
          'Bias Vmax')
xlab2 <- c("Absolute bias log(M)",
          "Absolute bias CR",
          'Absolute log(A50)',
          'Absolute log(beta)',
          'Absolute catch bias',
          'Absolute bias Vmax')

text.size <- 8 
Pnum <- 0
for (x in seq_along(funs)) {
  Phake <- funs[[x]]('Pacific_hake')
  Arrow <- funs[[x]]('Arrowtooth_Flounder')
  Swarehou <- funs[[x]]('Silver_warehou')
  POP <- funs[[x]]('Pacific_ocean_perch')
  
  DF <- bind_rows(Phake, Arrow, Swarehou, POP)
  DF$Species <- gsub("_", " ", DF$Species)
  DF$Species <- tolower(DF$Species)
  DF$Species <- paste0(toupper(substr(DF$Species, 1, 1)), substr(DF$Species, 2, nchar(DF$Species)))
  
  # ---- Relative Change ----
  Pnum <- Pnum + 1
  dat <- DF %>% filter(Years %in% c(Y1, Y2)) %>% 
    select(sim, bias, Years, Catch=stRelC, Depletion=relD, Species) %>% 
    tidyr::gather("key", "value", 4:5)
  dat$bias <- round(dat$bias, 2)
  
  Ylimits <- dat %>% group_by(Years, Species, key, bias) %>%
    summarize(upper=quantile(value, 0.95)) %>% ungroup() %>%
    summarize(max=max(upper))

  Pout <- ggplot(dat,
                 aes(x=as.factor(bias), y=value, fill=as.factor(Years))) +
    facet_grid(key~Species, scales = 'free') +
    geom_hline(yintercept = 1, linetype=3, col="darkgray") +
    geom_hline(yintercept = 0.5, linetype=4, col="darkgray") +
    geom_hline(yintercept = 2, linetype=4, col="darkgray") +
    geom_hline(yintercept = 0.25, linetype=5, col="darkgray") +
    geom_hline(yintercept = 4, linetype=5, col="darkgray") +
    geom_boxplot(outlier.shape = NA) + 
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(ylim=c(0, Ylimits$max)) +
    theme_classic() + theme(strip.background = element_blank()) + 
    theme(axis.title = element_text(size=text.size),
          axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90))+
    labs(fill="Year", x=xlab[x], y="Relative Change") +
    scale_fill_manual(values=cols1)
  
  filename <- paste0("Figures/SFigure", Pnum, ".png")
  ggsave(Pout, filename = filename, dpi=600, width=180, height=160, units="mm")
  
  # --- Expected Loss with over-estimation -----
  Pnum <- Pnum + 1
  
  dat <- DF %>% group_by(abs_bias, Species, Years, sim) %>% 
    summarise(Depletion=mean(relD[bias==bias[2]]/relD[bias==bias[1]]),
              Catch=mean(stRelC[bias==bias[2]]/stRelC[bias==bias[1]])) %>%
    filter(Years %in% c(Y1, Y2)) %>% tidyr::gather("key", "value", 5:6)
  
  Ylimits <- dat %>% group_by(Years, Species, key, abs_bias) %>%
    summarize(upper=quantile(value, 0.95)) %>% ungroup() %>%
    group_by(key) %>%
    summarize(max=max(upper))
  
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
    scale_y_continuous(expand=c(0,0),labels=scaleFUN) +
    coord_cartesian(ylim=c(0, min(Ylimits$max[1],10))) +
    theme_classic() + theme(strip.background = element_blank()) +
    theme(axis.title = element_text(size=text.size),
          axis.text = element_text(size=text.size-1),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text=element_text(size=text.size),
          strip.text.y=element_blank(),
          legend.position = 'top',
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    labs(fill="Year", x=xlab2[x], y="Catch Ratio") +
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
    scale_y_continuous(expand=c(0,0),labels=scaleFUN) +
    coord_cartesian(ylim=c(0, min(Ylimits$max[2],10))) +
    theme_classic() + theme(strip.background = element_blank()) +
    theme(axis.title = element_text(size=text.size),
          axis.text = element_text(size=text.size-1),
          axis.text.x = element_text(angle=90),
          strip.text.x = element_blank(),
          strip.text.y=element_blank(),
          strip.text=element_text(size=text.size),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    labs(fill="Year", x=xlab2[x], y="Depletion Ratio") +
    scale_fill_manual(values=cols2) + guides(fill=FALSE)
  
  leg <- cowplot::get_legend(P1)
  Plots <- cowplot::plot_grid(P1 + theme(legend.position="none"),
                              P2 + theme(legend.position="none"),
                              align = 'vh',
                              ncol=1, nrow=2,
                              hjust =-1)
  Pout <- cowplot::plot_grid(leg, Plots, rel_heights = c(0.2, 2), nrow=2)
  filename <- paste0("Figures/SFigure", Pnum, ".png")
  ggsave(Pout, filename = filename, dpi=600, width=180, height=160, units="mm")
  
}

