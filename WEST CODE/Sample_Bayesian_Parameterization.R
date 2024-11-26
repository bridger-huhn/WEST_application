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
library(progress)
rm(list = ls())




setwd("/Users/bridgerhuhn/Documents/Research/GYE_Endemics/")

### read in Data ####
dbConn <- dbConnect(SQLite(),"DATA/DataBase/EndemicPlantDataBase.db", overwrite = TRUE)

## list data tables
dbListTables(dbConn)


#read in phi2 data data
photo1 <- dbReadTable(dbConn,"PhotosynQ") %>% mutate(sppDateTrip = paste0(spp,"_",Date,Trip))
bpen <- dbReadTable(dbConn,"BlowPenComp")

#close db connection
dbDisconnect(dbConn)
###

###get bpen as comp for both bo sites

bpen <- rbind(bpen %>% mutate(Trip = "blowPen1"),bpen %>%
                mutate(Trip = "blowPen2")) %>% 
  mutate(sppDateTrip = paste0(spp,"_",Date,Trip))

photo1 <- rbind(bpen,photo1)
###
PhiCols <- names(photo1)[grep("Phi2",names(photo1))]
targetCol <- "Leaf.Temp.Differential"
photo1$plant.ID <- 1:(nrow(photo1))

### get c hab spp
temp <- photo1 %>% filter(hab == "c")

tempCHab <- photo1 %>%
  filter(sppDateTrip %in% unique(temp$sppDateTrip)) 

#### fit light response curves for each individual plant ####

# pivot data into long formate for easier regression fitting ####
photo <- tempCHab %>%
  select(spp,Date,rarity,hab,plant.ID,Trip,all_of(PhiCols),all_of(targetCol)) %>%
  pivot_longer(all_of(PhiCols)) %>% 
  separate(name, into = c("parameter", "PPFD"), sep = "_") %>% 
  mutate(PPFD = as.numeric(PPFD)) %>% 
  mutate(dateTrip = paste0(Date,Trip))
####


### set up a vector of error IDs
# Initialize vectors to store results and errors

params2 <- data.frame()
errors2 <- data.frame(plant.ID = character(), error.message = character(), stringsAsFactors = FALSE)

# Loop through each unique plant ID and assign curve to each individual plant measurement
for (ID in unique(photo$plant.ID)) {
  print(ID)
  # Filter data for the current plant ID and the specified parameter
  temp <- photo %>% filter(parameter == 'Phi2' & plant.ID == ID)
  
  # Try fitting the model, capturing any errors
  model <- try(nlsLM(
    value ~ ((a - k ) * exp(PPFD * (-1*b))) + k,
    data = temp,
    start = list(a = 0.7,  b = .0012, k = .08),
    lower = c(a = 0, b = 0, k = 0),
    upper = c(a = .9, b = .01, k = .35),
    control = nls.lm.control(maxiter = 1024)
  ), silent = TRUE)
  
  # Check if an error occurred
  if (inherits(model, "try-error")) {
    # If there was an error, store the plant ID and the error message in errors2
    error_message <- attr(model, "condition")$message
    errors2 <- rbind(errors2, data.frame(plant.ID = ID, error.message = error_message, stringsAsFactors = FALSE))
  } else {
    # If the model succeeded, extract the parameters and store them
    temp <- data.frame(plant.ID = ID, 
                       a = coef(model)[['a']],
                       b = coef(model)[['b']],
                       k = coef(model)[['k']])
    params2 <- rbind(params2, temp)
  }
}

### join the parameters and errors with their fluorescence measurements
photoParams <- left_join(params2,tempCHab ,by = c("plant.ID")) %>%
  select(Date, spp, plant.ID, hab, Trip,rarity, a,b,k) %>%
  mutate(dateTrip = paste0(Date,"_",Trip)) 


errors2 <- right_join(photo,errors2)


### set up bayesian model hyper parameters
iter = 8000
warmup = 1000
chains = 2
thin = 4

### this shows how many data points are in each Trip to check if there is enought to do analysis on
summaryPhoto <- photoParams %>% group_by(Date,Trip) %>% summarise(count =n())


#where to save models to 
fp <-"MODELS/ParamsModels/ByHabPhi_simple/"


#### A+B parameter analysis ####

for (trip in unique(photoParams$Trip)) {
  
  temp <- photoParams %>% filter(Trip == trip) 
  if(length(unique(temp$Date))==1){
    amod <- brm(bf(a ~ hab),
                data = temp,
                family = "beta",
                iter = iter,
                warmup = warmup,
                chains = chains,
                thin = thin,
                seed=1991,
                control = list(adapt_delta = .9)
                
    )
    
    print(paste0(trip,"---->alpha<---"))
    
    bmod <- brm(bf(b ~ hab),
                data = temp,
                family = "lognormal",
                iter = iter,
                warmup = warmup,
                chains = chains,
                thin = thin,
                seed=1991,
                control = list(adapt_delta = .9)
    )
    saveRDS(amod,paste0(fp,trip,"_a_beta_hab.RDS"))
    saveRDS(bmod,paste0(fp,trip,"_b_lognorm_hab.RDS"))
  }else{
    amod <- brm(bf(a ~ (1|Date)+hab),
                data = temp,
                family = "beta",
                iter = iter,
                warmup = warmup,
                thin = thin, 
                chains = chains,
                seed=1991,
                control = list(adapt_delta = .9)
    )
    print(paste0(trip,"---->alpha<---"))
    bmod <- brm(bf(b ~ (1|Date)+hab),
                data = temp,
                family = "lognormal",
                iter = iter,
                warmup = warmup,
                thin = thin,
                chains = chains,
                seed=1991,
                control = list(adapt_delta = .9)
    )
    saveRDS(amod,paste0(fp,trip,"_a_beta_hab.RDS"))
    saveRDS(bmod,paste0(fp,trip,"_b_lognorm_hab.RDS"))
  }
}


#### K parameters ####


for (trip in unique(photoParams$Trip)) {
  print(paste0("K_PARAMS------>",trip ,"<-----"))
  temp <- photoParams %>% filter(Trip == trip) 
  
  ## for single date
  if(length(unique(temp$Date))==1){
    kmod <- brm(bf(k ~ hab),
                data = temp,
                family = "zero_inflated_beta",
                iter = iter,
                warmup = warmup,
                chains = chains,
                thin = thin,
                seed=1991
                # prior = c(
                #   set_prior("normal(0.1, 0.05)", class = "Intercept", dpar = "mu")
                # ),
                # control = list(adapt_delta = .9)
    )
    
    saveRDS(kmod,paste0(fp,trip,"_k_ZIB_hab.RDS"))
  }else{
    kmod <- brm(bf(k ~ (1|Date)+hab),
                data = temp,
                family = "zero_inflated_beta",
                iter = iter,
                warmup = warmup,
                chains = chains,
                thin = thin,
                seed=1991
                # prior = c(
                #   set_prior("normal(0.1, 0.05)", class = "Intercept", dpar = "mu")
                # ),
                # control = list(adapt_delta = .9)
    )
    
    saveRDS(kmod,paste0(fp,trip,"_k_ZIB_hab.RDS"))
  }
  
  
  
}




