require(brms)
library(RSQLite)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(bayestestR)
require(tidybayes)
require(plotly)
require(minpack.lm)
require(loo)
require(cowplot)
library(egg)


rm(list = ls())



setwd("/Users/bridgerhuhn/Documents/Research/GYE_Endemics/")

##load in toSciNames function
source("CODE/functions/toSciNames.R") #### assigns scientific names to our data tags
source("CODE/functions/renameSite.R") ### this fixes sites to match the name of the rare plant of interest

fp <- "DATA/SummarizedModels/" ### file path for draws
fp2 <- "DATA/SummarizedHabModels/" ### file path for Habitat models
plotfp <- "FIGURES/PaperPlots/" ### where to save the plots to
fp_hab <- "Hab_PLOTS/" ### where to save the habitat plots to 


### read in data


## For Rarity models: read in all posterior draws and summaries and gets sci names
files <- paste0(fp,list.files(fp))
for (file in files) {
  
  assign(gsub(".RDS","", basename(file)),toSciNames(readRDS(file)))

}

## For Hab models: read in all draws and summaries and gets sci names
files <- paste0(fp2,list.files(fp2))
for (file in files) {
  
  assign(paste0("Hab_",gsub(".RDS","", basename(file)))
         ,toSciNames(readRDS(file)))
  
}

### read in data ##

### parameters
PhotoParams <- readRDS("DATA/fluorescenceParameters.RDS") %>% filter(Phi2_250>.1)



## set the credible interval to be plotted in all of the graphs


CI <- .19
Chi <- 1-(CI/2)
Clow <- 1-(Chi)


### set plot theme ###

cust_theme <- function(legend.pos = 'none'){
  theme(legend.position = legend.pos,
        legend.justification = c("right", "top"),
        panel.background = element_blank(), 
        panel.grid.major = element_line(color = "grey75"),  # Major grid lines in light grey
        panel.grid.minor = element_line(color = "grey83"),   # Remove minor grid lines
        #axis.line = element_line(color = "black"),     # Add axis lines
        axis.ticks = element_blank(),                  # Remove axis ticks
        axis.title = element_text(size = 12),          # Customize axis title
        plot.title = element_text(size = 14, hjust = 0.5), # Center the title
        axis.text = element_text(size = 10)  )}


### 1:1 plots #####
fp_1v1 <- "1v1_PLOTS/" ### where to save these plots

#### Nitrogen content ####
ggplot(NSimVDat_summaries, aes(x = sim_mean, y = dat_mean, color = rarity)) +
  geom_abline(slope = 1,intercept = 0)+
  geom_errorbar(aes(xmax = sim_qi.66, xmin = sim_qi.33), size = 1.5)+
  geom_errorbar(aes(xmax = sim_qi.90,xmin = sim_qi.10),size = .5)+
  geom_errorbar(aes(ymax = dat_qi.66, ymin = dat_qi.33), size = 1.5)+
  geom_errorbar(aes(ymax = dat_qi.90,ymin = dat_qi.10),size = .5)+
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black")+
  xlab(expression(Simulated~Leaf~Nitrogen~"%"))+
  ylab(expression(Measured~Leaf~Nitrogen~"%"))+
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  cust_theme(legend.pos = c(.2,.9))

## C13
ggplot(C13SimvDat, aes(x = sim_mean, y = dat_mean, color = rarity)) +
  geom_abline(slope = 1,intercept = 0)+
  geom_errorbar(aes(xmax = sim_qi.66, xmin = sim_qi.33), size = 1.5)+
  geom_errorbar(aes(xmax = sim_qi.90,xmin = sim_qi.10),size = .5)+
  geom_errorbar(aes(ymax = dat_qi.66, ymin = dat_qi.33), size = 1.5)+
  geom_errorbar(aes(ymax = dat_qi.90,ymin = dat_qi.10),size = .5)+
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black")+
  xlab(expression(Simulated~delta^13*C ~ "\U2030"))+
  ylab(expression(Measured~delta^13*C ~ "\U2030"))+
  xlim(-31,-21.5)+
  ylim(-31,-21.5)+
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41"))+
  cust_theme(legend.pos = c(.2,.9))




## Phi_PSII 

### fixes rarity variables 
PhiModVDat$rarity[which(PhiModVDat$rarity == "r")] <- "e"

Hab_Phi2DrawsVDat_summaries$hab[which(Hab_Phi2DrawsVDat_summaries$hab == "r")] <- "e"

NPQDrawsVDat_summaries$rarity[which(NPQDrawsVDat_summaries$rarity=="r")]<- "e"

Hab_NPQDrawsVDat_summaries$hab[which(Hab_NPQDrawsVDat_summaries$hab=="r")]<- "e"



Phi2_1v1_rarity <- ggplot(PhiModVDat, aes(x = draws_mean, y = dat_mean, shape = rarity,color = rarity, group=(interaction(Trip,rarity)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0)+
  xlab(expression(Simulated~Mean~~phi~PSII))+
  ylab(expression(Measured~Mean~~phi~PSII))+
  xlim(0,.65)+
  ylim(0,.65)+
  labs(color = "Species Type",shape = "Species Type")+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41"))+
  cust_theme(legend.pos = c(.4,.9))

Phi2_1v1_hab <- ggplot(Hab_Phi2DrawsVDat_summaries, aes(x = draws_mean, y = dat_mean, shape = hab,color = hab, group=(interaction(Trip,hab)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0)+
  xlab(expression(Simulated~Mean~~phi~PSII))+
  ylab(expression(Measured~Mean~~phi~PSII))+
  labs(color = "Species Type",shape = "Species Type")+
  xlim(0,.65)+
  ylim(0,.65)+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41")) +
  cust_theme(legend.pos = c(.4,.9))


NPQ_1v1_rarity<- ggplot(NPQDrawsVDat_summaries, aes(x = draws_mean, y = dat_mean, shape = rarity,color = rarity, group=(interaction(Trip,rarity)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0, alpha = .5)+
  xlab(expression(Simulated~Mean~~phi~NPQ))+
  ylab(expression(Measured~Mean~~phi~NPQ))+
  xlim(0,.9)+
  ylim(0,.9)+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41")) +
  labs(color = "Habitat Type",shape = "Habitat Type")+
  cust_theme(legend.pos = c(.4,.9))

NPQ_1v1_hab<- ggplot(Hab_NPQDrawsVDat_summaries, aes(x = draws_mean, y = dat_mean, shape = hab,color = hab, group=(interaction(Trip,hab)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0, alpha = .5)+
  xlab(expression(Simulated~Mean~~phi~NPQ))+
  ylab(expression(Measured~Mean~~phi~NPQ))+
  xlim(0,.9)+
  ylim(0,.9)+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41")) +
  labs(color = "Habitat Type",shape = "Habitat Type")+
  cust_theme(legend.pos = c(.4,.9))

Phi2_NPQ_1v1 <- plot_grid(plotlist = list(NPQ_1v1_hab,NPQ_1v1_rarity,
          Phi2_1v1_hab,Phi2_1v1_rarity))
ggsave(paste0(plotfp,fp_1v1, "Phi2_NPQ_1v1.pdf"),Phi2_NPQ_1v1, width = 3000, height = 3000, units = "px")
### dat V Simulation Plots####



# N_percent
ggplot(NSimVDat, aes(y = value, x = Trip, color = rarity, fill = simulation, shape = simulation)) +
  stat_pointinterval(side = "both",position = "dodgejust") +
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  xlab(expression(Leaf~Nitrogen~"%"))+
  ylab("Site")+
  guides(x =  guide_axis(angle = -45))+
  theme_minimal()

# Delta_C 
ggplot(C13_Draws, aes(y = value, x = Trip, color = rarity, fill = simulation, shape = simulation)) +
  ylim(-33, -20) +
  stat_pointinterval(side = "both",position = "dodgejust") +
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  xlab(expression(delta^13*C ~ "\U2030"))+
  ylab("Site")+
  guides(x =  guide_axis(angle = -45))+
  theme_minimal()

## Phi_PSII
ggplot(Phi_Draws, aes(x = parn, y = Phi, color = rarity, fill = points, shape = points))+
  stat_halfeye(position="dodge")+
  xlab(expression(PPFD ~ (mu*mol~m^-2~s^-1)))+
  ylab(expression(Measured~phi*PSII))+
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  ggtitle("Modeled vs Measured PSII Efficiency")+
  scale_shape_manual(values=c(2, 8))+
  theme_minimal()

### Phi_NPQ
ggplot(NPQDrawsVDat, aes(x = parn, y = NPQ, color = rarity, fill = points, shape = points))+
  stat_halfeye(position="dodge")+
  xlab(expression(PPFD ~ (mu*mol~m^-2~s^-1)))+
  ylab(expression(Measured~phi*NPQ))+
  scale_color_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  ggtitle("Modeled vs Measured\nNon-Phototochemical Quenching")+
  scale_shape_manual(values=c(2, 8))+
  cust_theme(legend.pos = "right")+
  ylim(0,1)


### a params

ggplot(NPQParamsDraws, aes(x = a , y = Trip, fill = rarity))+
  geom_violin()+
  scale_fill_manual(values = c("r" = "#5b76c6", "c" = "#f5ba41")) +
  geom_vline(xintercept = 0)+
  xlab(expression(pi[NPQ]))+
  ggtitle(paste0("Effect of Rarity on P_zi Estimate \n",CI,"% Credible Intervals"))+
  cust_theme(legend.pos = "bottom")



# Rarity Parameters ####
fp_rarity <- "Rarity_PLOTS/" 
fp_hab <- "Hab_PLOTS/"
#### NPQ ####

NPQ_rarityparams <- NPQ_rarityparams %>% filter(Trip != "spread" & Trip !="pilgrim" )


kRare <- ggplot(NPQ_rarityparams %>% filter(param == "k_rarityr"),
                aes(x = Estimate, y = Trip)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), size = 1, color = "#5b76c6") +  # Set color for error bars
  geom_vline(xintercept = 0) +
  ylab(NULL) +
  xlab(NULL)+
  ggtitle(bquote(kappa[Phi~NPQ]))+
  theme(axis.text.y = element_blank()) +
  cust_theme()


bRare <- ggplot(NPQ_rarityparams %>% filter(param == "b_rarityr"),
                aes(x = Estimate , y = Trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle(expression(beta[Phi~NPQ]))+
  theme(axis.text.y = element_blank())+
  cust_theme()

aRare <- ggplot(NPQ_rarityparams %>% filter(param == "a_rarityr"),
                aes(x = Estimate , y = Trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),
                size = 1,color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(expression(alpha[Phi~NPQ]))+
  cust_theme()

NPQ_rarity <- ggarrange(aRare, bRare,
                        ncol = 2,
                        nrow=1,
                        widths = c(1, 1),
                        top = bquote("Effect Of Endemism On The Mean of"))
#### Phi ####
Phi_rarityparams<- PhiII_rarity_params %>% filter(trip != "spread" & trip !="pilgrim" )

kRare <- ggplot(Phi_rarityparams %>% filter(param == "k_rarityr"), aes(x = Estimate , y = trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  xlab(expression(pi[Phi~PSII[endemic]]))+
  ylab(NULL)+
  ggtitle(expression(kappa[Phi~PSII]))+
  theme(axis.text.y = element_blank())+
  cust_theme()

kZiRare <- ggplot(Phi_rarityparams %>% filter(param == "k_zi_rarityr"),
                  aes(x = Estimate, color = trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  xlab(expression(pi[Phi~PSII[endemic]]))+
  ylab(NULL)+
  ggtitle(expression(kappa[Phi~PSII]))+
  theme(axis.text.y = element_blank())+
  cust_theme()



bRare <- ggplot(Phi_rarityparams %>% filter(param == "b_rarityr"),
                aes(x = Estimate , y = trip, color = trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1, color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  xlab(expression(mu[Phi~PSII[endemic]]))+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote(beta[Phi~PSII]))+
  theme(axis.text.y = element_blank())+
  cust_theme()

aRare <- ggplot(Phi_rarityparams %>% filter(param == "a_rarityr"),
                aes(x = Estimate , y = trip))+
  geom_point(color = "#5b76c6")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1, color = "#5b76c6")+
  geom_vline(xintercept = 0)+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(expression(alpha[Phi~PSII]))+
  cust_theme()

Phi2_rarity<- ggarrange(aRare, bRare, ncol = 2, nrow=1, widths = c(1, 1),
                        top = bquote("Effect Of Endemism On The Mean of"))

phiVNPQ <- plot_grid(NPQ_rarity, Phi2_rarity, nrow = 2)

ggsave(paste0(plotfp,fp_rarity,"PhiVNPQ_rarity.png"), phiVNPQ,width = 2900,height = 2900, units = "px")

# Habitat Plots ####

#1:1###

Phi2_1v1 <- ggplot(Hab_Phi2DrawsVDat_summaries, aes(x = draws_mean, y = dat_mean, shape = hab,color = hab, group=(interaction(Trip,hab)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0, alpha = .6)+
  xlab(expression(Simulated~Mean~~phi~PSII))+
  ylab(expression(Measured~Mean~~phi~PSII))+
  xlim(0,.8)+
  ylim(0,.8)+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41")) +
  cust_theme(legend.pos = c(.2,.9))

ggsave(Phi2_1v1,
       filename = paste0(plotfp,fp_1v1,"HabPhi2_1v1",".pdf"))

NPQ_1v1<- ggplot(Hab_NPQDrawsVDat_summaries, aes(x = draws_mean, y = dat_mean, shape = hab,color = hab, group=(interaction(Trip,hab)))) +
  geom_point(alpha = .5)+
  geom_line(alpha = .5)+
  geom_abline(slope = 1,intercept = 0, alpha = .6)+
  xlab(expression(Simulated~Mean~~phi~NPQ))+
  ylab(expression(Measured~Mean~~phi~NPQ))+
  xlim(0,.85)+
  ylim(0,.85)+
  scale_color_manual(values = c("e" = "#5b76c6", "c" = "#f5ba41")) +
  cust_theme(legend.pos = c(.2,.9))

ggsave(NPQ_1v1,
       filename = paste0(plotfp,fp_1v1,"HabNPQ_1v1",".pdf"))

ggarrange(NPQ_1v1, Phi2_1v1)
#### rarity Params ####
##### Phi ####
Hab_Phi2_habparams$Trip[which(Hab_Phi2_habparams$Trip == "spread")] <- "Spread Creek Site"

plot_Phi_a_hab <- ggplot(Hab_Phi2_habparams %>% filter(param == "a_habr"),
                         aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  geom_vline(xintercept = 0)+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle(bquote(alpha[Phi~PSII]))+
  cust_theme()

plot_Phi_b_hab <- ggplot(Hab_Phi2_habparams  %>% filter(param == "b_habr"),
                         aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  geom_vline(xintercept = 0)+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle(bquote(beta[Phi~PSII]))+
  theme(axis.text.y = element_blank())+
  cust_theme()

plot_Phi_k_hab <- ggplot(Hab_Phi2_habparams  %>% filter(param == "k_habr"),
                         aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  geom_vline(xintercept = 0)+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle(bquote(kappa[Phi~PSII]))+
  theme(axis.text.y = element_blank())+
  cust_theme()

### put plots into one plot
plot_Phi_rarity_hab <- ggarrange(plot_Phi_a_hab, plot_Phi_b_hab,  ncol = 2, nrow=1, widths = c(1, 1),
                                 top = bquote("Effect Of Habitat On The Mean of"))

Hab_NPQ_habparams$Trip[which(Hab_NPQ_habparams$Trip == "spread")] <- "Spread Creek Site"

NPQ_a_habr <- ggplot(Hab_NPQ_habparams %>% filter(param == "a_habr"), aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  ylab(NULL)+
  geom_vline(xintercept = 0)+
  xlab(NULL)+
  ggtitle(bquote(alpha[Phi~NPQ]))+
  cust_theme()


NPQ_b_habr <- ggplot(Hab_NPQ_habparams  %>% filter(param == "b_habr"), aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  geom_vline(xintercept = 0)+
  ylab(NULL)+
  theme(axis.text.y = element_blank())+
  xlab(NULL)+
  ggtitle(bquote(beta[Phi~NPQ]))+
  cust_theme()

NPQ_k_habr<- ggplot(Hab_NPQ_habparams  %>% filter(param == "k_habr"), aes(x = Estimate , y = Trip))+
  geom_point(color = "#f5ba41")+
  ylab(NULL)+
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi),size = 1,color = "#f5ba41")+
  geom_vline(xintercept = 0)+
  theme(axis.text.y = element_blank())+
  xlab(NULL)+
  ggtitle(bquote(kappa[Phi~NPQ]))+
  cust_theme()

plot_NPQ_rarity_Hab<- ggarrange(NPQ_a_habr,NPQ_b_habr, ncol = 2, nrow=1, widths = c(1, 1))

Hab_rarity <- plot_grid(plot_Phi_rarity_hab, plot_NPQ_rarity_Hab, nrow = 2)

ggsave(paste0(plotfp,"Hab_rarity_params",".png"),Hab_rarity,width = 2900,height = 2900, units = "px")


### PARAMETER and DATA plots #####
m250<-lm(Phi2_250~PhiNPQ_250, data = PhotoParams )
m500<- lm(Phi2_500~PhiNPQ_500, data = PhotoParams)
m1000<- lm(Phi2_1000~PhiNPQ_1000, data = PhotoParams)
m1500<- lm(Phi2_1500~PhiNPQ_1500, data = PhotoParams)
m2000<- lm(Phi2_2000~PhiNPQ_2000, data = PhotoParams)

p250 <- ggplot(PhotoParams,aes(PhiNPQ_250,Phi2_250, color = PhiNO_250))+
  geom_point()+
  xlim(0,1)+
  ylim(0,.7)+
  xlab("")+
  ylab(expression(phi~"PSII"))+
  ggtitle("250")+
  geom_abline(slope = coef(m250)[2], intercept = coef(m250)[1])+
  scale_color_gradient(low =  "blue", high = "red") +
  cust_theme()

p500<- ggplot(PhotoParams,aes(PhiNPQ_500,Phi2_500, color = PhiNO_500))+
  geom_point()+
  geom_abline(slope = coef(m250)[2], intercept = coef(m250)[1])+
  xlim(0,1)+
  ylim(0,.7)+
  ylab(NULL)+
  xlab("")+
  ggtitle("500")+
  theme(axis.text.y = element_blank()) +
  scale_color_gradient(low =  "blue", high = "red") +
  cust_theme()

p1000<- ggplot(PhotoParams,aes(PhiNPQ_1000,Phi2_1000, color = PhiNO_1000))+
  geom_point()+
  geom_abline(slope = coef(m250)[2], intercept = coef(m250)[1])+
  xlim(0,1)+
  ylim(0,.7)+
  ylab(NULL)+
  xlab(expression(phi~"NPQ"))+
  ggtitle("1000")+
  theme(axis.text.y = element_blank()) +
  scale_color_gradient(low =  "blue", high = "red") +
  cust_theme()

p1500<-ggplot(PhotoParams,aes(PhiNPQ_1500,Phi2_1500, color = PhiNO_1500))+
  geom_point()+
  xlim(0,1)+
  ylim(0,.7)+
  geom_abline(slope = coef(m250)[2], intercept = coef(m250)[1])+
  ylab(NULL)+
  xlab("")+
  ggtitle("1500")+
  theme(axis.text.y = element_blank()) +
  scale_color_gradient(low =  "blue", high = "red") +
  cust_theme()

p2000<-ggplot(PhotoParams,aes(PhiNPQ_2000,Phi2_2000, color = PhiNO_2000))+
  geom_point()+
  geom_abline(slope = coef(m250)[2], intercept = coef(m250)[1])+
  xlim(0,1)+
  ylim(0,.7)+
  ylab(NULL)+
  xlab("")+
  ggtitle("2000")+
  theme(axis.text.y = element_blank()) +
  scale_color_gradient(low =  "blue", high = "red") +
  labs(color = expression(Phi~NO))+
  cust_theme(legend.pos = "right")


NOdeviance <- ggarrange(p250,p500,p1000,p1500,p2000,
                        nrow =1,
                        widths = c(1.25,1,1,1,1.7),
          bottom = expression(Phi~PSII))

title <- text_grob(expression("PPFD (" * mu * "mol " * m^{-2} * s^{-1} * ")"))

NOdeviance<- annotate_figure(NOdeviance,
                top = title)

ggsave(paste0(plotfp, "NOdeviance",".png"), NOdeviance)



ggplot(PhotoParams,aes(b,PhiNO_250, color = PhiNPQ_250))+
  geom_point()+

  ggtitle(expression("2000 PPFD (" * mu * "mol " * m^{-2} * s^{-1} * ")"))+
  scale_color_gradient(low =  "blue", high = "red") +
  cust_theme()

NOcols <- names(PhotoParams)[grep("PhiNO",names(PhotoParams))]
t <-PhotoParams %>%
  select(plant.ID,all_of(NOcols)) %>%
  pivot_longer(cols = all_of(NOcols)) %>%
  separate(name, sep = "_",into = c("trash","PAR")) %>%
  select(-trash) %>% 
  mutate(PAR = as.numeric(PAR))


t<- left_join(t, PhotoParams, by = "plant.ID")
ggplot(t, aes(PAR, value,group = plant.ID, color = NPQ_b, alpha = .01))+geom_line()+theme_minimal()+
  scale_color_gradient(low =  "blue", high = "red")


temp <- t %>% group_by(plant.ID) %>%
  summarize(b = mean(b),
            NPQ_b = mean(NPQ_b),
            a = mean(a),
            NPQ_a = mean(NPQ_a),
            k = mean(k),
            NPQ_k = mean(NPQ_k),
            min = min(value),
            max = max(value),
            sd = sd(value)) %>%
  mutate(spread = max-min)

### HAB vs Endemism models ####

Hab_Phi2_habparams$Trip[which(Hab_Phi2_habparams$Trip == "spread")] <- "Spread Creek Site"
Phi_rarityparams<- PhiII_rarity_params %>% filter(trip != "pilgrim")
Phi_rarityparams[which(Phi_rarityparams$Trip == "spread")] <- "Spread Creek Site"

habVendParams <- Phi_rarityparams[Phi_rarityparams$trip %in% Hab_Phi2_habparams$Trip,]

Hab_Phi2_habparams <- Hab_Phi2_habparams %>%
  separate(param, into = c("param","trash"),sep = "_") %>% select(-trash)


habVendParams <- habVendParams %>%
  filter(param != "k_zi_rarityr") %>%
  separate(param, into = c("param","trash"),sep = "_") %>% select(-trash, -temp) %>% 
  rename(Trip = trip,
         phiEstimate = Estimate,
         phiCIlow = CIlow,
         phiCIhi = CIhi)

habVendParams<- left_join(Hab_Phi2_habparams, habVendParams, by = c("Trip", "param"))


a_habVend <- ggplot(habVendParams %>% filter(param == "a"),
                aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), linewidth = 1, color = "#f5ba41") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black", nudge_y = .004)+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-0.06,.15)+
  ylim(-0.06,.15)+
  ggtitle(bquote(alpha[Phi~PSII]))+
  ylab(bquote("Effect of Endemism on the Mean of"~alpha))+
  xlab(bquote("Effect of Endemic Habitat on the Mean of"~alpha))+
  cust_theme()

b_habVend <- ggplot(habVendParams %>% filter(param == "b"),
                    aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), size = 1, color = "#f5ba41") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black", nudge_y = .00004)+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  # xlim(-0.0012,.001)+
  # ylim(-0.0012,.001)+
  ggtitle(bquote(beta[Phi~PSII]))+
  ylab(bquote("Effect of Endemism on the Mean of"~beta))+
  xlab(bquote("Effect of Endemic Habitat on the Mean of"~beta))+
  cust_theme()

k_habVend <- ggplot(habVendParams %>% filter(param == "k"),
                    aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), size = 1, color = "#5b76c6") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black")+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle(bquote(kappa[Phi~PSII]))+
  cust_theme()

EndVHab<- ggarrange(a_habVend, b_habVend,  ncol = 2, nrow=1, widths = c(1, 1))
ggsave("FIGURES/PaperPlots/Hab_V_Endemism.pdf",EndVHab, width = 3000, height = 1500, units = "px")

### NPQ 

Hab_NPQ_habparams$Trip[which(Hab_NPQ_habparams$Trip == "spread")] <- "Spread Creek Site"
NPQ_rarityparams<- NPQ_rarityparams %>% filter(Trip != "pilgrim")
NPQ_rarityparams$Trip[which(NPQ_rarityparams$Trip == "spread")] <- "Spread Creek Site"

habVendParams <- NPQ_rarityparams[NPQ_rarityparams$Trip %in% Hab_NPQ_habparams$Trip,]

Hab_NPQ_habparams <- Hab_NPQ_habparams %>%
  separate(param, into = c("param","trash"),sep = "_") %>% select(-trash)


habVendParams <- habVendParams %>%
  filter(param != "k_zi_rarityr") %>%
  separate(param, into = c("param","trash"),sep = "_") %>% select(-trash, -trash) %>% 
  rename(Trip = Trip,
         phiEstimate = Estimate,
         phiCIlow = CIlow,
         phiCIhi = CIhi)

habVendParams<- left_join(Hab_NPQ_habparams, habVendParams, by = c("Trip", "param"))


a_habVend_NPQ <- ggplot(habVendParams %>% filter(param == "a"),
                    aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), linewidth = 1, color = "#8CC63F") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black", nudge_y = .004)+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-0.16,.15)+
  ylim(-0.16,.15)+
  ggtitle(bquote(alpha[Phi~NPQ]))+
  ylab(bquote("Effect of Endemism on the Mean of"~alpha))+
  xlab(bquote("Effect of Endemic Habitat on the Mean of"~alpha))+
  cust_theme()

b_habVend_NPQ <- ggplot(habVendParams %>% filter(param == "b"),
                    aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), size = 1, color = "#8CC63F") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black", nudge_y = .00004)+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-0.0015,.001)+
  ylim(-0.0015,.001)+
  ggtitle(bquote(beta[Phi~NPQ]))+
  ylab(bquote("Effect of Endemism on the Mean of"~beta))+
  xlab(bquote("Effect of Endemic Habitat on the Mean of"~beta))+
  cust_theme()

k_habVend <- ggplot(habVendParams %>% filter(param == "k"),
                    aes(x = Estimate, y = phiEstimate)) +
  geom_point(color = "#5b76c6") +  # Set color for points
  geom_errorbar(aes(xmin = CIlow, xmax = CIhi), size = 1, color = "#5b76c6") + 
  geom_errorbar(aes(ymin = phiCIlow, ymax = phiCIhi), size = 1, color = "#5b76c6") + 
  geom_label(aes(label = Trip), size = 2.5,label.size = .02, color = "black")+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle(bquote(kappa[Phi~PSII]))+
  cust_theme()

EndVHab<- ggarrange(a_habVend, b_habVend,
                    a_habVend_NPQ, b_habVend_NPQ,
                    ncol = 2, nrow=2, widths = c(1, 1))
ggsave("FIGURES/PaperPlots/Hab_V_Endemism.pdf",EndVHab, width = 3000, height = 3000, units = "px")

