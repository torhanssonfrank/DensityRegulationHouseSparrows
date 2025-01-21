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
growth_factors<-read.table("Data/growth_factors.txt", sep = ";")

#Set the recruit year as 0.
adult_fitness<- adult_fitness %>% 
  mutate(Least_age = Least_age-1) %>% 
  mutate(Least_age_sq = Least_age^2)

x<-adult_fitness %>% 
  dplyr::select(Location, in_out) %>% 
  unique() %>% 
  summarise(sum(in_out))

stopifnot(x==6)

#r2S-models####

##meta-population N_m model ####
adult_fitness$w<-adult_fitness$r2s

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd")) #You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

#fit a model without the inout interaction with density.
bernt_zNR_nb_brm_no_in_out_slope<-brm(w ~ zNR + in_out + Least_age_sq+ Least_age  * sex + (1|ID) + (1+zNR|Location)
                                      + (1|Year_Location),adult_fitness,iter = 25000,thin = 30,
                                      cores = 4, negbinomial(link = "log", link_shape = "log"))

saveRDS(bernt_zNR_nb_brm_no_in_out_slope, file("Workspace backup/bernt_zNR_nb_brm_no_in_out_slope.rds"))

#fit main model
bernt_zNR_nb_brm_noindslope<-brm(w ~ zNR*in_out + Least_age_sq+ Least_age  * sex + (1|ID) + (1+zNR|Location)
                                 + (1|Year_Location),adult_fitness,iter = 25000,thin = 30,
                                 cores = 4, negbinomial(link = "log", link_shape = "log"))

saveRDS(bernt_zNR_nb_brm_noindslope, file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))

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
      ind_isl_mod <- brm(w ~ gamma +Least_age_sq+ Least_age * sex
                         + (1|ID) 
                         + (1|Year),ind_isl_data, iter = ni,thin = nt,chains = nc,
                         cores = nc, prior = bs_prior, negbinomial(link = "log", link_shape = "log"))
    } else if(distribution == "bernoulli"){
      ind_isl_mod <- brm(w ~ gamma + Least_age_sq+Least_age  * sex
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

adult_fitness$gamma<-adult_fitness$zNR
distribution = "negbin"
ni = 12000
nt = 30
nc = 4

results_single_islands<-single_island_analysis(adult_fitness,distribution,ni,nt,nc,bs_prior)


saveRDS(results_single_islands, file("Workspace backup/bernt_zNR_single_island_mods.rds"))

##log N model ####

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd")) #You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

#fit a model without the inout interaction with density.
bernt_n_nb_brm_no_in_out_slope<-brm(w ~ n+in_out + Least_age_sq + Least_age * sex + (1|ID) + (1+n|Location)
                                    + (1|Year_Location),adult_fitness,iter = 12000,thin = 30,
                                    cores = 4, prior = bs_prior ,negbinomial(link = "log", link_shape = "log"))
saveRDS(bernt_n_nb_brm_no_in_out_slope, file("Workspace backup/bernt_n_nb_brm_no_in_out_slope.rds"))


#main model
bernt_n_nb_brm_noindslope<-brm(w ~ n*in_out + Least_age_sq + Least_age * sex + (1|ID) + (1+n|Location)
                               + (1|Year_Location),adult_fitness,iter = 12000,thin = 30,
                               cores = 4, prior = bs_prior ,negbinomial(link = "log", link_shape = "log"))
saveRDS(bernt_n_nb_brm_noindslope, file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))

#Recruit and nestling production models ####
##recruit N_m ####

adult_fitness$w<-adult_fitness$N_Recruits

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd"))#You cannot put lb (lower bound) or ub(upper bound) on SD's in BRMS. But the lower bound is 0 by default.

rec_zNR_nb_brm_noindslope<-brm(w ~ zNR*in_out + Least_age_sq + Least_age * sex+ (1|ID) + (1+zNR|Location)
                               + (1|Year_Location), adult_fitness,iter = 25000,thin = 30, cores = 4,prior=bs_prior, negbinomial(link = "log", link_shape = "log"))

saveRDS(rec_zNR_nb_brm_noindslope, file("Workspace backup/rec_zNR_nb_brm_noindslope.rds"))

##nestling production N_m ####

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd"))


nestlings_zNR_nb_brm_noindslope<-brm(N_nestlings ~ zNR*in_out + Least_age_sq + Least_age * sex + (1|ID) + (1+zNR|Location)
                                     + (1|Year_Location),nestling_producers,iter = 12000,thin = 30,
                                     cores = 4, prior = bs_prior, negbinomial(link = "log", link_shape = "log"))
saveRDS(nestlings_zNR_nb_brm_noindslope, file("Workspace backup/nestlings_zNR_nb_brm_noindslope.rds"))

#Sensitivities####

#remove NA rows
growth_factors<-na.omit(growth_factors)


#remove rows with 0
growth_factors<-growth_factors %>% 
  filter(!if_any(c(8:11),~ . < 1))


growth_factors<-growth_factors %>% 
  mutate(across(8:11, log))

growth_factors<-growth_factors %>% 
  filter(lambda<4)

growth_factors<-growth_factors %>% 
  mutate(n = log(N))

bs_prior <- c(set_prior("normal(0,2)", class = "b"),
              set_prior("normal(0,1)", class = "sd"))


growth_mod_lambda_brm<-brm(lambda ~ n + Adult_survival + 
                             Fledgling_to_rec + 
                             Nestling_to_fledge + 
                             Nestlings_born +
                             (1|Year) + (n|Location), 
                           data = growth_factors,
                           family = Gamma(link = "log"),
                           prior = bs_prior,
                           iter = 60000,
                           thin = 60,
                           cores = 4)
saveRDS(growth_mod_lambda_brm, "Workspace backup/growth_mod_lambda_brm.rds")

growth_mod_lambda_brm_inout<-brm(lambda ~ n + 
                                   in_out * Adult_survival + 
                                   in_out * Fledgling_to_rec + 
                                   in_out * Nestling_to_fledge + 
                                   in_out * Nestlings_born +
                                   (1|Year) + (n|Location), 
                                 data = growth_factors,
                                 family = Gamma(link = "log"),
                                 prior = bs_prior,
                                 iter = 60000,
                                 thin = 60,
                                 cores = 4)
saveRDS(growth_mod_lambda_brm_inout, "Workspace backup/growth_mod_lambda_brm_inout.rds")


max_ages<- nestling_producers %>% 
  group_by(ID) %>% 
  mutate(Least_age = Least_age + 1) %>% 
  mutate(max_age = max(Least_age)) %>% 
  mutate(lrs = sum(N_nestlings)) %>% 
  mutate(nestlings_per_year = lrs/max_age) %>% 
  dplyr::select(ID, Location, max_age, nestlings_per_year, lrs, sex) %>% 
  unique()

descdist(max_ages$max_age, boot = 500, discrete = T)

age_mod_brm<-brm(max_age~nestlings_per_year + (1|Location),
                 data = max_ages,
                 family = poisson,
                 prior = bs_prior,
                 iter = 20000,
                 thin = 30,
                 cores = 4)

saveRDS(age_mod_brm, "Workspace backup/max_age_yearly_nestlingprod_brm.rds")
