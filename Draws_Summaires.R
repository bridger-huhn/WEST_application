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

rm(list = ls())


setwd("/Users/bridgerhuhn/Documents/Research/GYE_Endemics/")


### set the Credible intervals you wish to see on the graphs ###
CI <- .19
Chi <- 1-(CI/2)
Clow <- 1-(Chi)



### a_ models ####

fp <- ("MODELS/ParamsModels/NPQHabBest/")

files <- list.files(fp)[grep("_a_", list.files(fp))]





### this extracts draws form the models
for (file in files) {
  print(file)
  mod <- readRDS(paste0(fp,file))
  
  ### make a dataframe of the summary
  temp <- data.frame(posterior_summary(mod, probs = c(Clow,Chi)))
  
  ###
  temp <- temp %>%
    mutate(param = rownames(temp)) %>%
    filter(param != "lp__" & param != "lprior") %>%
    filter(param != "Intercept" & param != "Intercept_phi" & param != "Intercept_zi") %>% 
    mutate(Trip = gsub("_a_beta_hab.RDS","",file)) 
  
  
  
  common <- temp %>% filter(param == "b_Intercept"|param == "b_zi_Intercept" | param == "b_phi_Intercept") %>% 
    select(1,3,4,5)
  rare <- temp %>% filter(param == 'b_habr'|param == "b_zi_habr" | param == " b_phi_hab") %>% 
    select(1,3,4,5)
  rare <- data.frame(Estimate = inv_logit_scaled(rare$Estimate[1]+common$Estimate[1]),
                     CIlow =inv_logit_scaled(rare[1,2]+common[1,1]),
                     CIhi = inv_logit_scaled(rare[1,3]+common[1,1]),
                     param = rare$param) 
  
  common_params <- common$param
  common <- common %>%
    select(-param) %>% mutate(across(everything(),~inv_logit_scaled(.x))) %>% 
    mutate(param = common_params)
  temp <- data.frame(Estimate = rare$Estimate[1]-common$Estimate[1],
                     CIlow = rare[1,2]-common[1,1],
                     CIhi = rare[1,3]-common[1,1],
                     param = rare$param) %>% 
    mutate(trip = gsub("NPQ_","",gsub("_a_beta_hab.RDS","",file)))
  
  ### we need different simulations because the some models don't have date as a random efffet
  dates <- mod[["basis"]][["levels"]][["Date"]]
  if(is.null(dates)){
    ## extract posterior samples as predictions
    
    
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = data.frame(hab = factor(c("r","c"))),
                                                         ndraws = 3500))
    names(posterior_samples) <- c("r","c")
    
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab")%>% 
      mutate(Date = "SingleDate")
    
    
    
  }else{
    ### made data with multiple dates
    newdata = data.frame(hab = factor(c(rep("r",length(dates)),rep("c",length(dates)))),
                         Date = factor(rep(dates,2)))
    post_names <- paste0(newdata$hab,"_",newdata$Date)
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = newdata,
                                                         ndraws = 3500))
    
    #fix names of posterior samples
    names(posterior_samples) <- post_names
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab") %>% 
      separate(hab, into = c("hab", "Date"), sep = "_") 
    
  }
  
  posterior_samples$Trip <- gsub("NPQ_","",gsub("_a_beta_hab.RDS","",file))
  
  
  if(file == files[1]){
    a_post_draws <- posterior_samples
    a_habparams <- temp
  }else{
    a_post_draws <- rbind(a_post_draws,posterior_samples)
    a_habparams <- rbind(a_habparams, temp)
  }
}

#### summarize a_ Draws ####
a_simDat_Summaries <- a_post_draws %>% group_by(hab, Trip) %>% 
  summarize(mean = mean(value),
            CIhi = quantile(value, probs = Chi),
            CIlow = quantile(value, probs = Clow))




### b model draws ####


files <- list.files(fp)[grep("_b_", list.files(fp))]


for (file in files) {
  print(file)
  mod <- readRDS(paste0(fp,file))
  temp <- data.frame(posterior_summary(mod, probs = c(Clow,Chi)))
  temp <- temp %>%
    mutate(param = rownames(temp)) %>%
    filter(param == "b_habr" | param == "b_Intercept") %>%
    mutate(Trip = gsub("NPQ_","",gsub("_b_lognorm_hab.RDS","",file)))
  common <- temp %>% filter(param == "b_Intercept") %>% 
    select(1,3,4)
  rare <- temp %>% filter(param == 'b_habr') %>% 
    select(1,3,4)
  rare <- data.frame(Estimate = inv_logit_scaled(rare$Estimate[1]+common$Estimate[1]),
                     CIlow =inv_logit_scaled(rare[1,2]+common[1,1]),
                     CIhi = inv_logit_scaled(rare[1,3]+common[1,1])) 
  common <- common %>% mutate(across(everything(),~inv_logit_scaled(.x)))
  temp <- data.frame(Estimate = rare$Estimate[1]-common$Estimate[1],
                     CIlow = rare[1,2]-common[1,2],
                     CIhi = rare[1,3]-common[1,3]) %>% 
    mutate(trip = gsub("NPQ_","",gsub("_b_lognorm_hab.RDS","",file)))
  
  ### we need different simulations because the some models don't have date as a random efffet
  dates <- mod[["basis"]][["levels"]][["Date"]]
  if(is.null(dates)){
    ## extract posterior samples as predictions
    
    
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = data.frame(hab = factor(c("r","c"))),
                                                         ndraws = 3500) )
    names(posterior_samples) <- c("r","c")
    
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab") %>% 
      mutate(Date = "SingleDate")
    
    
    
  }else{
    ### made data with multiple dates
    newdata = data.frame(hab = factor(c(rep("r",length(dates)),rep("c",length(dates)))),
                         Date = factor(rep(dates,2)))
    post_names <- paste0(newdata$hab,"_",newdata$Date)
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = newdata,
                                                         ndraws = 3500))
    
    #fix names of posterior samples
    names(posterior_samples) <- post_names
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab") %>% 
      separate(hab, into = c("hab", "Date"), sep = "_") 
    
  }
  
  posterior_samples$Trip <- gsub("NPQ_","",gsub("_b_lognorm_hab.RDS","",file))
  
  
  if(file == files[1]){
    b_post_draws <- posterior_samples
    b_habparams <- temp
  }else{
    b_post_draws <- rbind(b_post_draws,posterior_samples)
    b_habparams <- rbind(b_habparams, temp)
  }
}

### k model draws ####


files <- list.files(fp)[grep("_k_", list.files(fp))]


for (file in files) {
  print(file)
  mod <- readRDS(paste0(fp,file))
  temp <- data.frame(posterior_summary(mod, probs = c(Clow,Chi)))
  temp <- temp %>%
    mutate(param = rownames(temp)) %>%
    filter(param == "b_habr" | param == "b_Intercept"|
             param == "b_zi_habr"| param =="b_zi_Intercept") %>%
    mutate(Trip = gsub("NPQ_","",gsub("_k_ZIB_hab.RDS","",file))) 
  common <- temp %>% filter(param =="b_Intercept") %>% 
    select(1,3,4,5)
  rare <- temp %>% filter(param =="b_habr") %>% 
    select(1,3,4,5)
  rare <- data.frame(Estimate = inv_logit_scaled(rare$Estimate+common$Estimate),
                     CIlow =inv_logit_scaled(rare[,2]+common[,1]),
                     CIhi = inv_logit_scaled(rare[,3]+common[,1]),
                     param = rare$param) 
  common[,-4]<-inv_logit_scaled(common[,-4])
  temp <- data.frame(Estimate = rare$Estimate-common$Estimate,
                     CIlow = rare[,2]-common$Estimate,
                     CIhi = rare[,3]-common$Estimate,
                     Trip = gsub("NPQ_","",gsub("_k_ZIB_hab.RDS","",file)),
                     param = rare$param) 
  
  ### we need different simulations because the some models don't have date as a random efffet
  dates <- mod[["basis"]][["levels"]][["Date"]]
  if(is.null(dates)){
    ## extract posterior samples as predictions
    
    
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = data.frame(hab = factor(c("r","c"))),
                                                         ndraws = 3500) )
    names(posterior_samples) <- c("r","c")
    
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab") %>% 
      mutate(Date = "SingleDate")
    
    
    
  }else{
    ### made data with multiple dates
    newdata = data.frame(hab = factor(c(rep("r",length(dates)),rep("c",length(dates)))),
                         Date = factor(rep(dates,2)))
    post_names <- paste0(newdata$hab,"_",newdata$Date)
    posterior_samples <- as.data.frame(posterior_predict(mod,
                                                         newdata = newdata,
                                                         ndraws = 3500))
    
    #fix names of posterior samples
    names(posterior_samples) <- post_names
    
    posterior_samples <- posterior_samples %>%
      pivot_longer(everything(), names_to = "hab") %>% 
      separate(hab, into = c("hab", "Date"), sep = "_") 
    
  }
  
  posterior_samples$Trip <- gsub("NPQ_","",gsub("_k_ZIB_hab.RDS","",file))
  
  
  if(file == files[1]){
    k_post_draws <- posterior_samples
    k_habparams <- temp
  }else{
    k_post_draws <- rbind(k_post_draws,posterior_samples)
    k_habparams <- rbind(k_habparams, temp)
  }
}


### if zi varies with rarity
# for (file in files) {
#   print(file)
#   mod <- readRDS(paste0(fp,file))
#   temp <- data.frame(posterior_summary(mod, probs = c(Clow,Chi)))
#   temp <- temp %>%
#     mutate(param = rownames(temp)) %>%
#     filter(param == "b_rarityr" | param == "b_Intercept"|
#              param == "b_zi_rarityr"| param =="b_zi_Intercept") %>%
#     mutate(Trip = gsub("NPQ_","",gsub("_k_ZIB_rarity.RDS","",file))) 
#   common <- temp %>% filter(param =="b_Intercept") %>% 
#     select(1,3,4,5)
#   rare <- temp %>% filter(param =="b_rarityr") %>% 
#     select(1,3,4,5)
#   rare <- data.frame(Estimate = inv_logit_scaled(rare$Estimate+common$Estimate),
#                      CIlow =inv_logit_scaled(rare[,2]+common[,1]),
#                      CIhi = inv_logit_scaled(rare[,3]+common[,1]),
#                      param = rare$param) 
#   common[,-4]<-inv_logit_scaled(common[,-4])
#   temp <- data.frame(Estimate = rare$Estimate-common$Estimate,
#                      CIlow = rare[,2]-common$Estimate,
#                      CIhi = rare[,3]-common$Estimate,
#                      Trip = gsub("NPQ_","",gsub("_k_ZIB_rarity.RDS","",file)),
#                      param = rare$param) 
#   
#   if(file == files[1]){
#     k_b_rarityparams <- temp
#   }else{
#     k_b_rarityparams <- rbind(k_b_rarityparams, temp)
#   }
# }

### join Data ####
habparams <- rbind(k_habparams %>% mutate(param = "k_habr"),
                      #k_b_habparams %>% mutate(param = "k_b_habr"),
                      a_habparams%>%
                        mutate(param = "a_habr")  %>% rename(Trip = trip),
                      b_habparams %>% mutate(param = "b_habr") %>% 
                        rename(Trip = trip))

habparams$Trip[grep("yermo_north",habparams$Trip)]<-"yermoNorth"
habparams$Trip[grep("yermo_south",habparams$Trip)]<-"yermoSouth"
habparams$Trip<- sub("_.*", "", habparams$Trip)



### join posteriors
# test to see if hab and Trip columns match before 
which(a_post_draws$hab != b_post_draws$hab) ### should == 0
which(k_post_draws$hab != b_post_draws$hab) ### should == 0
which(k_post_draws$Trip != b_post_draws$Trip) ### should == 0
which(a_post_draws$Trip != b_post_draws$Trip) ### should == 0
which(a_post_draws$Date != b_post_draws$Date) ### should == 0
which(k_post_draws$Date != b_post_draws$Date) ### should == 0

draws <- cbind(a_post_draws %>% select(-Date, -hab,-Trip) %>% rename(a = value),
               b_post_draws %>% rename(b = value) %>% select(-Date,-hab, -Trip),
               k_post_draws %>% rename(k=value))

rm(list = ls()[grep("temp", ls())])
rm(list = ls()[grep("mod", ls())])


### bind data together

##make a function for NPQ
NPQcurve <- function(a,b,k, ppfd){
  ppfd <- (a-(1-k))*exp(ppfd*-b)+(1-k)
  return(ppfd)
}

## calculate NPQ at different light levels

draws <- draws %>% 
  mutate(NPQ_250 = NPQcurve(a,b,k,250)) %>% 
  mutate(NPQ_500 = NPQcurve(a,b,k,500)) %>% 
  mutate(NPQ_1000 = NPQcurve(a,b,k,1000)) %>% 
  mutate(NPQ_1500 = NPQcurve(a,b,k,1500)) %>% 
  mutate(NPQ_2000 = NPQcurve(a,b,k,2000))

NPQCols <- names(draws)[grep("NPQ", names(draws))]
NPQParamDraws <- draws
tdraws <- draws %>% select(hab,Date,Trip,NPQCols) %>%
  mutate(plantID = 1:nrow(draws)) %>%
  pivot_longer(all_of(NPQCols), names_to = "PAR", values_to = "NPQ") %>%
  separate(PAR,sep = "_", into = c("trash","PAR")) %>%
  select(-trash) %>%
  mutate(parn = as.numeric(PAR)) %>%
  mutate(PAR = fct_reorder(PAR, parn)) %>%
  mutate(points = "Modeled")




### read in photosynQ data #####

dbConn <- dbConnect(SQLite(),"DATA/DataBase/EndemicPlantDataBase.db")

photo <- dbReadTable(dbConn,"photosynQ")  %>% mutate(sppDateTrip = paste0(spp,"_",Date,Trip))

### get common hab spp
temp <- photo %>% filter(hab == "c")

photo <- photo %>%
  filter(sppDateTrip %in% unique(temp$sppDateTrip)) 




NPQNames <- names(photo)[grep("NPQ_",names(photo))]

photo <- photo %>%
  mutate(plantID = (1:nrow(photo))) %>% 
  select(plantID, hab,Trip,Date,all_of(NPQNames)) %>%
  pivot_longer(cols = all_of(NPQNames), names_to = "PAR", values_to = "NPQ") %>%
  separate(PAR, sep = "_", into = c("trash","PAR")) %>% 
  select(-trash) %>% mutate(parn = as.numeric(PAR)) %>% 
  mutate(PAR = fct_reorder(PAR, parn)) %>% 
  mutate(points = "Data")



### phi data summary
photo_dat_summary <- photo %>%  
  group_by(Trip, hab,PAR) %>% 
  summarize(mean = mean(NPQ, na.rm = TRUE),
            qi.10 = quantile(NPQ, probs = .1, na.rm = TRUE),
            qi.33 = quantile(NPQ, probs = .33, na.rm = TRUE),
            qi.66 = quantile(NPQ,probs = .66, na.rm = TRUE),
            qi.90 = quantile(NPQ, probs = .9, na.rm = TRUE)) %>% 
  mutate(data = "data") 
## put a tag on names of data
names(photo_dat_summary)[-c(1:3,ncol(photo_dat_summary))]<- paste0("dat_",names(photo_dat_summary)[-c(1:3,ncol(photo_dat_summary))])

### summarize draws
photo_draws <- draws %>% select(hab,Date,Trip,all_of(NPQCols)) %>% 
  mutate(plantID = 1:nrow(draws)) %>% 
  pivot_longer(NPQCols, names_to = "PAR", values_to = "NPQ") %>% 
  separate(PAR, into = c("trash","PAR")) %>% 
  select(-trash) %>% 
  mutate(parn = as.numeric(PAR)) %>% 
  mutate(PAR = fct_reorder(PAR, parn)) %>% 
  mutate(points = "Modeled")

draws_datVmodeled <- rbind(photo,photo_draws)

### summarize draws
photo_draws_summary <- photo_draws %>% 
  group_by(Trip,hab, PAR) %>% 
  summarize(mean = mean(NPQ),
            qi.10 = quantile(NPQ, probs = .1),
            qi.33 = quantile(NPQ, probs = .33),
            qi.66 = quantile(NPQ,probs = .66),
            qi.90 = quantile(NPQ, probs = .9)) %>% 
  mutate(data = "sim")


## put a tag on names of data
names(photo_draws_summary)[-c(1:3,ncol(photo_draws_summary))]<- paste0("draws_",names(photo_draws_summary)[-c(1:3,ncol(photo_draws_summary))])


### bind together
photoModVDat <- merge(photo_draws_summary, 
                      photo_dat_summary,
                      by = c("Trip","hab","PAR"))



### rename data sets for saving them
NPQDrawsVDat <- draws_datVmodeled
NPQDrawsVDat_summaries <- photoModVDat
NPQ_habparams <- habparams
### save data sets
### save draws and compiled data sets ####
fp <- "DATA/SummarizedHabModels/"
saveRDS(NPQDrawsVDat, paste0(fp,"NPQDrawsVDat",".RDS"))
saveRDS(NPQDrawsVDat_summaries, paste0(fp,"NPQDrawsVDat_summaries",".RDS"))
saveRDS(NPQ_habparams, paste0(fp,"NPQ_habparams",".RDS"))
saveRDS(NPQParamDraws, paste0(fp,"NPQParamsDraws.RDS"))

