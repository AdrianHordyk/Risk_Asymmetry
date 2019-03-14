
OM <- readRDS(file.path("OMs", paste0(SpeciesList[2], ".rdata")))

OM@h <- c(0.7, 0.7)
msys <- CalcMSYs(OM)

YieldC <- CalcYieldCurve(OM) %>% mutate(Yield=Yield/msys$Yield)

Mbias <- exp(c(-0.5, 0.5))
msys_2 <- CalcMSYs(OM, M=OM@M * Mbias[1])
msys_3 <- CalcMSYs(OM, M=OM@M * Mbias[2])

PlotList <- MakePlots(OM, msys, msys_2, msys_3)
df <- PlotList$df
df2 <- PlotList$df2

theme <- theme(axis.title=element_text(size=14), axis.text=element_text(size=12))

width <- 4.5; height <- 4
ggsave("Figures/Presentation/Y_F1.png", PlotList$Y_F + theme, width=width, height=height)
ggsave("Figures/Presentation/Y_D1.png", PlotList$Y_D + theme, width=width, height=height)
ggsave("Figures/Presentation/F_D1.png", PlotList$F_D + theme, width=width, height=height)


# Add MSY points
P1 <- PlotList$Y_F + 
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=Y), size=5) +
  geom_line(data=data.frame(x=c(msys$F, msys$F), y=c(-Inf, 1)), 
            aes(x=x, y=y), linetype=2) + theme

P2 <- PlotList$Y_D + 
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=B, y=Y), size=5) +
  geom_line(data=data.frame(x=c(msys$SB_SB0, msys$SB_SB0), y=c(-Inf, 1)), 
            aes(x=x, y=y), linetype=2) + theme


P3 <- PlotList$F_D  + theme + 
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=B), size=5) +
  geom_line(data=data.frame(x=c(msys$F, msys$F), y=c(-Inf, msys$SB_SB0)), 
            aes(x=x, y=y), linetype=2) +
  geom_line(data=data.frame(x=c(-Inf, msys$F), y=c(msys$SB_SB0, msys$SB_SB0)), 
            aes(x=x, y=y), linetype=2)

ggsave("Figures/Presentation/Y_F2.png", P1, width=width, height=height)
ggsave("Figures/Presentation/Y_D2.png", P2, width=width, height=height)
ggsave("Figures/Presentation/F_D2.png", P3, width=width, height=height)


# Add biases
size <- 2

OM2 <- OM; OM2@M <- OM@M * exp(-0.5)
YieldC2 <- CalcYieldCurve(OM2) %>% mutate(Yield=Yield/max(Yield))

OM3 <- OM; OM3@M <- OM@M * exp(0.5)
YieldC3 <- CalcYieldCurve(OM3) %>% mutate(Yield=Yield/max(Yield))

P1 <- PlotList$Y_F +
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=Y), size=5) + theme +
  scale_x_continuous(limits = c(0, 0.45))
ggsave("Figures/Presentation/Y_F3a.png", P1, width=width, height=height) 


P1 <- PlotList$Y_F +
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=Y), size=5) + 
  geom_line(data=YieldC2, aes(x=F, y=Yield), size=lsize, linetype=2)  + theme +
  geom_vline(data=msys_2, aes(xintercept = F), linetype=3) +
  scale_x_continuous(limits = c(0, 0.45)) 

ggsave("Figures/Presentation/Y_F3b.png", P1, width=width, height=height)


P1 <- PlotList$Y_F +
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=Y), size=5) + 
  geom_line(data=YieldC3, aes(x=F, y=Yield), size=lsize, linetype=3) +
  geom_line(data=YieldC2, aes(x=F, y=Yield), size=lsize, linetype=2) +
  geom_vline(data=msys_3, aes(xintercept = F), linetype=3) +
  theme +
  scale_x_continuous(limits = c(0, 0.45))
ggsave("Figures/Presentation/Y_F3c.png", P1, width=width, height=height)


P1 <- PlotList$F_D +
  geom_point(data=filter(df, bias=="Unbiased"), aes(x=F, y=B), size=5) + 
  geom_point(data=filter(df, bias=="Negative"), aes(x=F, y=B), size=5, shape=15) + 
  geom_point(data=filter(df, bias=="Positive"), aes(x=F, y=B), size=5, shape=17) +
  theme +
  scale_x_continuous(limits = c(0, 0.45))
ggsave("Figures/Presentation/F_D3.png", P1, width=width, height=height)


width <- 5.5; height <- 4
df$bias <- factor(df$bias, levels=c("Unbiased", "Negative", "Positive"), ordered = TRUE)
P1 <- PlotList$Y_F + 
  geom_point(data=df, aes(x=F, y=Y, shape=bias), size=5) + theme 
ggsave("Figures/Presentation/Y_F4.png", P1, width=width, height=height)

P2 <- PlotList$Y_D + 
  geom_point(data=df, aes(x=B, y=Y, shape=bias), size=5) + theme
ggsave("Figures/Presentation/Y_D4.png", P2, width=width, height=height)



lsize <- 1.25
PDF <- PlotList$PDF

P1 <- ggplot(filter(PDF, Year>50, bias=="Unbiased"), aes(x=Year, y=Biomass, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + 
  labs(x="Year", y="Spawning Depletion", shape="Bias", linetype='Bias') +
  expand_limits(y=c(0,1)) + theme + 
  geom_hline(yintercept = msys$SB_SB0, linetype=2) + 
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_1.png", P1, width=width, height=height)

P1 <- ggplot(filter(PDF, Year>50, bias!="Positive"), aes(x=Year, y=Biomass, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + 
  labs(x="Year", y="Spawning Depletion", shape="Bias", linetype='Bias') +
  expand_limits(y=c(0,1)) + theme +
  geom_hline(yintercept = msys$SB_SB0, linetype=2) + 
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_2.png", P1, width=width, height=height)

P1 <- ggplot(filter(PDF, Year>50), aes(x=Year, y=Biomass, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + theme +
  labs(x="Year", y="Spawning Depletion", shape="Bias", linetype='Bias') +
  expand_limits(y=c(0,1)) +
  geom_hline(yintercept = msys$SB_SB0, linetype=2) + 
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_3.png", P1, width=width, height=height)



P1 <- ggplot(filter(PDF, Year>50, bias=="Unbiased"), aes(x=Year, y=RelC, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + 
  labs(x="Year", y="Yield", shape="Bias", linetype='Bias') +
  geom_hline(yintercept = 1, linetype=2) + theme +
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_1a.png", P1, width=width, height=height)

P1 <- ggplot(filter(PDF, Year>50, bias!="Positive"), aes(x=Year, y=RelC, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + 
  labs(x="Year", y="Yield", shape="Bias", linetype='Bias') +
  geom_hline(yintercept = 1, linetype=2) +  theme +
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_2a.png", P1, width=width, height=height)

P1 <- ggplot(filter(PDF, Year>50), aes(x=Year, y=RelC, linetype=fct_inorder(bias))) + geom_line(size=lsize) +
  expand_limits(y=0,1) + theme_classic() + 
  labs(x="Year", y="Yield", shape="Bias", linetype='Bias') +
  geom_hline(yintercept = 1, linetype=2) +  theme +
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  guides(linetype=FALSE)
ggsave("Figures/Presentation/Proj_3a.png", P1, width=width, height=height)

# eg 20% bias + vs 20% -ve bias 


PDF %>% head()
temp <- PDF %>% filter(bias !="Unbiased")
temp$B <- temp$Biomass[temp$bias=="Positive"]/temp$Biomass[temp$bias=="Negative"]
temp$C <- temp$RelC[temp$bias=="Positive"]/temp$RelC[temp$bias=="Negative"]


P1 <- ggplot(filter(temp, Year>50), aes(x=Year, y=B)) + 
  geom_line(size=lsize) + expand_limits(y=c(0,1)) + 
  theme_classic() + theme +
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  labs(y='Depletion Ratio (postive/negative bias') +
  geom_hline(yintercept = 1, linetype=2)

ggsave("Figures/Presentation/Proj_4.png", P1, width=width, height=height)

P1 <- ggplot(filter(temp, Year>50), aes(x=Year, y=C)) + 
  geom_line(size=lsize) + expand_limits(y=c(0,1)) + 
  theme_classic() + theme +
  scale_x_continuous(labels=seq(0, 50, by=10), breaks=seq(50, 100, by=10)) +
  labs(y='Catch Ratio (postive/negative bias') +
  geom_hline(yintercept = 1, linetype=2)
ggsave("Figures/Presentation/Proj_4a.png", P1, width=width, height=height)

temp %>% select(Year, B)
temp %>% select(Year, C)


PComb <- grid_arrange_shared_legend(list(Y_F, Y_D, P_C, P_B), position = 'right') 
