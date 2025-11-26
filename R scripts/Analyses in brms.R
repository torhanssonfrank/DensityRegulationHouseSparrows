#BRMS models


library(arm)
library(rstan)
library(parallel)
library(tidyverse)
library(shinystan)
library(brms)
library(writexl)

adult_fitness<-read.table("Data/cleaned_density_dependence_data.txt", sep=";", header=TRUE)
nestling_producers <- read.table("Data/nestling_producers_2007_2014.txt", sep = ";")
demo <-read.table("Data/estimated_pop_sizes_nestling_to_adult_model_with_lurøy_onøy.txt", sep =";", header = T, encoding = "UTF-8")

demo_new_sub<-demo %>% 
  dplyr::select(Location, Year, N_corr) %>% 
  mutate(N_corr = round(N_corr))

adult_fitness<-adult_fitness %>% 
  left_join(demo_new_sub)

adult_fitness %>% 
  filter(is.na(N_corr))

island_N_corr_means<-adult_fitness %>% 
  dplyr::select(Year, Location, N_corr) %>% 
  unique() %>% 
  group_by(Location) %>% 
  summarise(Mean_N_corr_isl = mean(N_corr)) %>% 
  ungroup()

#Add the means back to the dataset
adult_fitness<- adult_fitness %>% 
  left_join(island_N_corr_means, by ="Location")

island_N_corr_means %>% 
  arrange(desc(Mean_N_corr_isl))

adult_fitness<-adult_fitness %>% 
  mutate(mc_N_corr = N_corr - Mean_N_corr_isl)

scaler <- 40 

adult_fitness<- adult_fitness %>% 
  mutate(zNR_new = mc_N_corr/scaler)

adult_fitness$n_new <- log(adult_fitness$N_corr)

nestling_producers<-nestling_producers %>% 
  left_join(demo_new_sub)


island_N_corr_means2<-nestling_producers %>% 
  dplyr::select(Year, Location, N_corr) %>% 
  unique() %>% 
  group_by(Location) %>% 
  summarise(Mean_N_corr_isl = mean(N_corr)) %>% 
  ungroup()

#Add the means back to the dataset
nestling_producers<- nestling_producers %>% 
  left_join(island_N_corr_means2, by ="Location")

island_N_corr_means2 %>% 
  arrange(desc(Mean_N_corr_isl))

nestling_producers<-nestling_producers %>% 
  mutate(mc_N_corr = N_corr - Mean_N_corr_isl)

scaler <- 40 

nestling_producers<- nestling_producers %>% 
  mutate(zNR_new = mc_N_corr/scaler)

#Set the recruit year as 0.
adult_fitness<-adult_fitness %>% 
  mutate(Least_age_mc = Least_age-mean(Least_age)) %>% 
  mutate(Least_age_mc_sq = Least_age_mc^2)

x<-adult_fitness %>% 
  dplyr::select(Location, in_out) %>% 
  unique() %>% 
  summarise(sum(in_out))

stopifnot(x==6)

nestling_producers<-nestling_producers %>% 
  mutate(Least_age_mc = Least_age-mean(Least_age)) %>% 
  mutate(Least_age_mc_sq = Least_age_mc^2)



#r2S-models####

##meta-population N_m model ####
adult_fitness$w<-adult_fitness$r2s

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd")) #You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

#fit a model without the inout interaction with density.
r2s_zNR_new_no_inout<-brm(w ~ zNR_new + in_out + Least_age_mc+ Least_age_mc_sq  * sex + (1|ID) + (1+zNR_new|Location)
                                      + (1|Year_Location),adult_fitness,iter = 8000,thin = 2,
                                      cores = 4, negbinomial(link = "log", link_shape = "log"))

saveRDS(r2s_zNR_new_no_inout, file("Workspace backup/r2s_zNR_new_no_inout.rds"))

#fit main model
r2s_zNR_new<-brm(w ~ zNR_new*in_out + Least_age_mc+ Least_age_mc_sq  * sex + (1|ID) + (1+zNR_new|Location)
                                 + (1|Year_Location),adult_fitness,iter = 8000,thin = 2,
                                 cores = 4, negbinomial(link = "log", link_shape = "log"))

saveRDS(r2s_zNR_new, file("Workspace backup/r2s_zNR_new.rds"))

##single island models####

single_island_analysis<-function(adult_fitness,distribution,ni,nt,nc,bs_prior){
  
  
  full_results <- list()
  
  ind_islands<-adult_fitness %>% 
    dplyr::select(Location) %>% 
    unique()
  brm
  for(i in 1:nrow(ind_islands)){
    ind_isl_data <- subset(adult_fitness, Location == ind_islands$Location[i])
    
    if(distribution == "negbin"){
      ind_isl_mod <- brm(w ~ gamma +Least_age_mc+ Least_age_mc_sq * sex
                         + (1|ID) 
                         + (1|Year),ind_isl_data, iter = ni,thin = nt,chains = nc,
                         cores = nc, prior = bs_prior, negbinomial(link = "log", link_shape = "log"))
    } else if(distribution == "bernoulli"){
      ind_isl_mod <- brm(w ~ gamma + Least_age_mc+Least_age_mc_sq  * sex
                         + (1|ID) 
                         + (1|Year),ind_isl_data, iter = ni,thin = nc,chains =nc,
                         cores = nc, prior = bs_prior, bernoulli(link = "logit"))
    } else{print("Error, distribution must be negbin or bernoulli!")}
    full_results[[i]] <- ind_isl_mod
    names(full_results)[i] <- unique(ind_isl_data$Location)
  }
  full_results
}

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd")) #You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

adult_fitness$gamma<-adult_fitness$zNR_new
distribution = "negbin"
ni = 6000
nt = 2
nc = 4

results_single_islands<-single_island_analysis(adult_fitness,distribution,ni,nt,nc,bs_prior)


saveRDS(results_single_islands, file("Workspace backup/r2s_zNR_single_island_mods_NEW.rds"))

##log N model ####

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd")) #You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

r2s_n_new<-brm(w ~ n_new*in_out + Least_age_mc + Least_age_mc_sq * sex + (1|ID) + (1+n_new|Location)
                               + (1|Year_Location),adult_fitness,iter = 4000,thin = 2,
                               cores = 4, prior = bs_prior ,negbinomial(link = "log", link_shape = "log"))
saveRDS(r2s_n_new, file("Workspace backup/r2s_n_new.rds"))

#Recruit and nestling production models ####
##recruit N_m ####

adult_fitness$w<-adult_fitness$N_Recruits

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd"))#You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

N_recruits_zNR_new<-brm(w ~ zNR_new*in_out + Least_age_mc + Least_age_mc_sq * sex+ (1|ID) + (1+zNR_new|Location)
                               + (1|Year_Location), adult_fitness,iter = 4000,thin = 2, cores = 4,prior=bs_prior, negbinomial(link = "log", link_shape = "log"))

saveRDS(N_recruits_zNR_new, file("Workspace backup/N_recruits_zNR_new.rds"))

##nestling production N_m ####

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd"))


nestlings_zNR_nb_brm_noindslope_new<-brm(N_nestlings ~ zNR_new*in_out + Least_age_mc_sq + Least_age_mc * sex + (1|ID) + (1+zNR_new|Location)
                                         + (1|Year_Location),nestling_producers,iter = 6000,thin = 6,
                                         cores = 4, prior = bs_prior, negbinomial(link = "log", link_shape = "log"))
saveRDS(nestlings_zNR_nb_brm_noindslope_new, file("Workspace backup/nestlings_zNR_nb_brm_noindslope_new.rds"))

