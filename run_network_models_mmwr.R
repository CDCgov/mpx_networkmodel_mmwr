# Here, we model the transmission of monkeypox through a population of 
# Gay, Bisexual, And other Men who have sex with Men using a 
# Stochastic Exponential Random Graph Approach. 
# This model estimates the impact of observed changes in partner acquisition behavior
# On monkeypox attack rates, and estimates the proportion of transmission to come
# from one time, casual, or main partnerships. 

# Written by Emily Pollock and Patrick Clay
# Last updated 8/23/2022

#load in libraries, set seed
library(here)
library(EpiModel)
library(EpiModelHIV)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(lhs)
library(kableExtra)
library(GGally)
library(tidyr)
set.seed(12345)

##############################

### Run Stergm ###

# load network structure
load(here("mpx_network_object.rda"))
nets <- est
rm(est)

# source modules for stergm
source(here("network_modules.R"))

#Set parameters
reps <- 60 #number of runs per parameter combination
cores <- detectCores() # don't do this on the cluster 
if (cores < reps) {ncores <- 6} #number of cores to use, can increase if desired
if (cores >= reps) {ncores <- reps} #number of cores to use


###########################################
#Run Low Transmission with Behavior Change
###########################################

# Set parameters that can change between stergm runs
freeparam1 <- 0.6 #probability of infection
freeparam2 <- 0.6 #probability of seeking treatment
freeparam3 <- 27  #average length of infectious period (days)
behavior.change.switch <- TRUE #Do individuals reduce number of one time partners?

# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 540, nsims=reps, ncores=ncores)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(nets, params, inits, controls)
#save results from this parameter set
save(sim,file="mmwr_behavechange_lowtrans.Rdata")
#remove all simulations with stochastic extinction
#defined as final # infected < initial # infected * 2
sim1<-sim #placeholder
for(i in 1:length(sim$epi)){
  placeholder <- sim$epi[[i]]  
  sim$epi[[i]] <- placeholder[which(sim1$epi$cuml.infs[365,] > 20)]
}
#extract medians
run_results_median <- extract_med_iqr(sim)
#save parameters
run_results_median$inf.prob <- freeparam1
run_results_median$treatment.prob <- freeparam2
run_results_median$days.til.recovery <- freeparam3
#save medians
write.csv(run_results_median, "mmwr_behavechange_lowtrans_noextinct_medians.csv")
#extract means
run_results_mean <- extract_mean_iqr(sim)
#save parameters
run_results_mean$inf.prob <- freeparam1
run_results_mean$treatment.prob <- freeparam2
run_results_mean$days.til.recovery <- freeparam3
#save means
write.csv(run_results_mean, "mmwr_behavechange_lowtrans_noextinct_means.csv")

##############################################
#Run Low Transmission without Behavior Change
##############################################

# Set parameters that can change between stergm runs
freeparam1 <- 0.6 #probability of infection
freeparam2 <- 0.6 #probability of seeking treatment
freeparam3 <- 27  #average length of infectious period (days)
behavior.change.switch <- FALSE #Do individuals reduce number of one time partners?

# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 540, nsims=reps, ncores=ncores)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(nets, params, inits, controls)
#save results from this parameter set
save(sim,file="mmwr_nobehavechange_lowtrans.Rdata")
#remove all simulations with stochastic extinction
#defined as final # infected < initial # infected * 2
sim1<-sim #placeholder
for(i in 1:length(sim$epi)){
  placeholder <- sim$epi[[i]]  
  sim$epi[[i]] <- placeholder[which(sim1$epi$cuml.infs[365,] > 20)]
}
#extract medians
run_results_median <- extract_med_iqr(sim)
#save parameters
run_results_median$inf.prob <- freeparam1
run_results_median$treatment.prob <- freeparam2
run_results_median$days.til.recovery <- freeparam3
#save medians
write.csv(run_results_median, "mmwr_nobehavechange_lowtrans_noextinct_medians.csv")
#extract means
run_results_mean <- extract_mean_iqr(sim)
#save parameters
run_results_mean$inf.prob <- freeparam1
run_results_mean$treatment.prob <- freeparam2
run_results_mean$days.til.recovery <- freeparam3
#save means
write.csv(run_results_mean, "mmwr_nobehavechange_lowtrans_noextinct_means.csv")


###########################################
#Run high Transmission with Behavior Change
###########################################

# Set parameters that can change between stergm runs
freeparam1 <- 0.9 #probability of infection
freeparam2 <- 0.6 #probability of seeking treatment
freeparam3 <- 27  #average length of infectious period (days)
behavior.change.switch <- TRUE #Do individuals reduce number of one time partners?

# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 540, nsims=reps, ncores=ncores)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(nets, params, inits, controls)
#save results from this parameter set
save(sim,file="mmwr_behavechange_hightrans.Rdata")
#remove all simulations with stochastic extinction
#defined as final # infected < initial # infected * 2
sim1<-sim #placeholder
for(i in 1:length(sim$epi)){
  placeholder <- sim$epi[[i]]  
  sim$epi[[i]] <- placeholder[which(sim1$epi$cuml.infs[365,] > 20)]
}
#extract medians
run_results_median <- extract_med_iqr(sim)
#save parameters
run_results_median$inf.prob <- freeparam1
run_results_median$treatment.prob <- freeparam2
run_results_median$days.til.recovery <- freeparam3
#save medians
write.csv(run_results_median, "mmwr_behavechange_hightrans_noextinct_medians.csv")
#extract means
run_results_mean <- extract_mean_iqr(sim)
#save parameters
run_results_mean$inf.prob <- freeparam1
run_results_mean$treatment.prob <- freeparam2
run_results_mean$days.til.recovery <- freeparam3
#save means
write.csv(run_results_mean, "mmwr_behavechange_hightrans_noextinct_means.csv")

##############################################
#Run high Transmission without Behavior Change
##############################################

# Set parameters that can change between stergm runs
freeparam1 <- 0.9 #probability of infection
freeparam2 <- 0.6 #probability of seeking treatment
freeparam3 <- 27  #average length of infectious period (days)
behavior.change.switch <- FALSE #Do individuals reduce number of one time partners?

# load in all other parameters
params <- param_msm()
# sequence of modules to run model
controls <- control_msm(nsteps = 540, nsims=reps, ncores=ncores)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(nets, params, inits, controls)
#save results from this parameter set
save(sim,file="mmwr_nobehavechange_hightrans.Rdata")
#remove all simulations with stochastic extinction
#defined as final # infected < initial # infected * 2
sim1<-sim #placeholder
for(i in 1:length(sim$epi)){
  placeholder <- sim$epi[[i]]  
  sim$epi[[i]] <- placeholder[which(sim1$epi$cuml.infs[365,] > 20)]
}
#extract medians
run_results_median <- extract_med_iqr(sim)
#save parameters
run_results_median$inf.prob <- freeparam1
run_results_median$treatment.prob <- freeparam2
run_results_median$days.til.recovery <- freeparam3
#save medians
write.csv(run_results_median, "mmwr_nobehavechange_hightrans_noextinct_medians.csv")
#extract means
run_results_mean <- extract_mean_iqr(sim)
#save parameters
run_results_mean$inf.prob <- freeparam1
run_results_mean$treatment.prob <- freeparam2
run_results_mean$days.til.recovery <- freeparam3
#save means
write.csv(run_results_mean, "mmwr_nobehavechange_hightrans_noextinct_means.csv")












