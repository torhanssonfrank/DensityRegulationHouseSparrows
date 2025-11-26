
#Survival from nestling to fledgeling and from fledgeling to recruit

library(Matrix)
library(tidyverse)
library(rstan)


data_observed<-read.table("Data/cleaned_density_dependence_data.txt", sep=";", header=TRUE)
demo <-read.table("Data/estimated_pop_sizes_nestling_to_adult_model_with_lurøy_onøy.txt", sep =";", header = T, encoding = "UTF-8")
pres_stage<-read.table("Data/juvenile_to_ad_surival_histories_long.txt", sep=";", header = T)

demo<-demo %>% 
  mutate(estimated_pop_size = round(N_corr))

nest_fledge<-pres_stage %>%
  dplyr::select(ID, Year, Location, Island, lifestage) %>% 
  unique() %>% 
  pivot_wider(names_from = lifestage, values_from = lifestage) %>% 
  relocate(nest, .before = fledge)

nest_fledge %>% 
  filter(Year == 2019 & (!is.na(recruit)|!is.na(adult))) #individuals born in 2019 still have a full CH.

CH_juv<-nest_fledge %>% 
  mutate(across(5:8, ~as.numeric(ifelse(is.na(.x), "0", "1")))) %>% 
  dplyr::select(-Year)
nrow(CH_juv)


#We just need the hatchyear and the born island

head(pres_stage)

year_juv<-pres_stage %>% 
  dplyr::select(ID, Year) %>% 
  unique()

year_juv %>% 
  group_by(ID) %>% 
  filter(n()>1)

nrow(year_juv)

island_juv<-pres_stage %>% 
  filter(lifestage == "nest") %>% 
  dplyr::select(ID, Island) %>%
  unique()
nrow(island_juv)


year_isl_juv<-island_juv %>% 
  left_join(year_juv) %>% 
  unite("FY", Island:Year)
nrow(year_isl_juv)



#Scale with roughly one sd (rounded to 40 inds)
scaler = 40

head(demo)
head(pres_stage)
##OBS! the island means should be calculated with the years used! ####

pop_year_island<-demo %>% 
  unite("FY",c(Island,Year), remove=F) 

pop_year_island<-pop_year_island %>% 
  filter(FY %in% year_isl_juv$FY)

pop_year_island %>% 
  group_by(Location) %>% 
  summarise(min(Year), max(Year))

pop_year_island<-pop_year_island %>% 
  dplyr::select(FY,Island, Year, estimated_pop_size) %>% 
  filter(!is.na(estimated_pop_size)) %>% 
  group_by(Island) %>% 
  mutate(mean_N_isl = mean(estimated_pop_size)) %>% 
  ungroup() %>% 
  mutate(zNR = (estimated_pop_size-mean_N_isl)/scaler) %>% 
  dplyr::select(FY, zNR) %>% 
  unique()


N<-year_isl_juv %>% 
  left_join(pop_year_island) %>% 
  dplyr::select(-FY)


Fl<-CH_juv %>%
  dplyr::select(ID, nest, fledge, recruit, adult) %>% 
  mutate(nest = 1, fledge = 0, recruit=0, adult = 0 )

R<-CH_juv %>%
  dplyr::select(ID, nest, fledge, recruit, adult) %>% 
  mutate(nest = 0, fledge = 1, recruit=0, adult = 0 )
A<-R %>% 
  mutate(nest = 0, fledge = 0, recruit=1, adult = 0 )

CH_juv<-as.data.frame(CH_juv)
A <- as.data.frame(A)

identical(year_juv[,1], island_juv[,1])
identical(CH_juv[,1], year_isl_juv[,1])
identical(CH_juv[,1], year_juv[,1])
identical(year_isl_juv[,1], island_juv[,1])
identical(A[,1],CH_juv[,1])
identical(N[,1],CH_juv[,1])

N <-N$zNR

Fl <- Fl[,2:ncol(Fl)]
R <- R[,2:ncol(R)]
A <- A[,2:ncol(A)]


Fl_p <- R
R_p <- A
A_p <- A
A_p[] <- 0
A_p$adult <- 1

#convert values to numerical starting at 1

island_juv$codeisl <- as.numeric(droplevels(as.factor(island_juv$Island)))

isl_name<-data_observed %>% 
  dplyr::select(Island, Location) %>% 
  unique()

name_code<-island_juv %>% 
  left_join(isl_name)

island<-as.numeric(droplevels(as.factor(island_juv$Island)))

FY <- as.numeric(droplevels(as.factor(year_isl_juv$FY)))

CH <- CH_juv[,4:ncol(CH_juv)]
nrow(CH_juv)
nrow(CH)

nFY <- max(FY)
in_out_short<-data_observed %>% 
  dplyr::select(Island, in_out) %>% 
  unique()

in_out<-island_juv %>% 
  left_join(in_out_short)

in_out <- in_out$in_out
length(in_out)


stan_data <- list(y =CH ,Fl = Fl, R =R, A=A,Fl_p = Fl_p, R_p = R_p, A_p = A_p,
                  in_out = in_out,  nind = dim(CH)[1], n_occasions =
                    dim(CH)[2], flok=island, nflok=max(island),  
                  N=N, FY=FY,nFY=nFY 
)


#Parameters monitored


params <- c("mu_p","mu_rec_p", "mu_phi","mu_rec_phi", "gamma","gamma_rec",
            "phi_fledge_interac_inout","phi_rec_interac_inout", 
            "gamma_fledge_interac_inout", "gamma_rec_interac_inout",
            "Sigma2_isl_phi_fledge" ,"Sigma2_isl_phi_rec",
            "Sigma2_isl_g_fledge","Sigma2_isl_g_rec",
            "Sigma2_YI_phi_fledge","Sigma2_YI_phi_rec",
            "Sigma2_YI_p_fledge", "Sigma2_YI_p_rec", 
            "cov_fledge","cov_rec",
            "isl_phi_fledge", "isl_phi_rec",
            "isl_g_fledge", "isl_g_rec")
## MCMC settings
ni <- 3000
nt <- 2
nc <- 4

mod <- stan_model("Stan scripts/CMR juvenile to adult one ind per isl full model.stan")


## Call Stan from R 
CMR_juvenile_survival_zNR_newpop <- sampling(mod,
                         data = stan_data, pars = params,
                         chains = nc, iter = ni, thin = nt,
                         cores=nc)

saveRDS(CMR_juvenile_survival_zNR_newpop, file("Workspace backup/CMR_juvenile_survival_zNR_newpop.rds"))
stan_res<-readRDS("Workspace backup/CMR_juvenile_survival_zNR_newpop.rds")


##Extract and print recruit gammas ####
uniq_islcode<-name_code %>% 
  dplyr::select(codeisl, Island, Location) %>% 
  unique() %>% 
  mutate(in_out = ifelse(Location %in% c("myken", "træna", "selvær", "lovund", "sleneset"), 0, 1)) %>% 
  arrange(codeisl)


rec_isl_gamma<-stan_res %>% 
  spread_draws(gamma_rec,gamma_rec_interac_inout, epsilonF_rec[variable, codeisl])


rec_isl_gamma<-rec_isl_gamma %>%
  filter(variable ==2)

head(rec_isl_gamma)

rec_isl_gamma %>% 
  median_qi(gamma_rec)

rec_isl_gamma %>% 
  median_qi(gamma_rec_interac_inout)

rec_isl_gamma2<-rec_isl_gamma %>% 
  left_join(uniq_islcode)

rec_isl_gamma2<-rec_isl_gamma2 %>% 
  mutate(gamma_fledge_to_rec_isl = ifelse(in_out == 1, gamma_rec + gamma_rec_interac_inout + epsilonF_rec, gamma_rec + epsilonF_rec)) %>% 
  group_by(Location) %>% 
  median_qi(gamma_fledge_to_rec_isl, .width = 0.9) %>% 
  ungroup()


saveRDS(rec_isl_gamma2, "Workspace backup/rec_isl_gamma.rds")


##Extract and print fledgling gammas ####
fledge_isl_gamma<-stan_res %>% 
  spread_draws(gamma,gamma_fledge_interac_inout, epsilonF_fledge[variable, codeisl])



fledge_isl_gamma<-fledge_isl_gamma %>%
  filter(variable ==2) %>% 
  ungroup()

head(fledge_isl_gamma)

fledge_isl_gamma %>% 
  median_qi(gamma)

fledge_isl_gamma %>% 
  median_qi(gamma_fledge_interac_inout)

fledge_isl_gamma2<-fledge_isl_gamma %>% 
  left_join(uniq_islcode)

fledge_isl_gamma2<-fledge_isl_gamma2 %>% 
  mutate(gamma_nest_to_fledge_isl = ifelse(in_out == 1, gamma + gamma_fledge_interac_inout + epsilonF_fledge, gamma + epsilonF_fledge)) %>% 
  group_by(Location) %>% 
  median_qi(gamma_nest_to_fledge_isl, .width = 0.9) %>% 
  ungroup()



saveRDS(fledge_isl_gamma2, "Workspace backup/fledge_isl_gamma.rds")





