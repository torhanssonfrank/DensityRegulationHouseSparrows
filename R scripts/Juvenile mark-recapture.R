
#Survival from nestling to fledgeling and from fledgeling to recruit

library(Matrix)
library(tidyverse)
library(rstan)


data_observed<-read.table("Data/cleaned_density_dependence_data.txt", sep=";", header=TRUE)
demo <-read.table("Data/demography.txt", sep =";", header = T)
pres_stage<-read.table("Data/juvenile_to_ad_surival_histories_long.txt", sep=";", header = T)


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



demo$estimated_pop_size[demo$Island == 22 &
                          demo$Year == 2014] <- 1 #Myken

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
  group_by(Islandname) %>% 
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
ni <- 6000
nt <- 30
nc <- 4


## Call Stan from R 
cjs_temp_raneff2 <- stan("Stan scripts/CMR juvenile to adult one ind per isl full model.stan",
                         data = stan_data, pars = params,
                         chains = nc, iter = ni, thin = nt,
                         cores=4)

saveRDS(cjs_temp_raneff2, file("Workspace backup/CMR_juvenile_survival_zNR.rds"))