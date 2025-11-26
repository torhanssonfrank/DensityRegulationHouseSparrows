#Survival from nestling to adult

library(Matrix)
library(tidyverse)
library(brms)
library(lme4)
library(rstan)
library(arm)
library(tidybayes)

pres<-read.delim("Data/presence_data_1994_2022.txt", sep =";", stringsAsFactors = F, encoding = "UTF-8")
sex_gen_corr<-read.table("Data/sex_genetically_corrected.txt", sep=";", header=TRUE)


##Create a stage column ####
pres_stage<-pres %>% 
  mutate(lifestage = ifelse(stage == "nest", "nest",
                            ifelse(stage != "nest" & Least_age == 0, "fledge",
                                   ifelse(Least_age == 1, "recruit",
                                          ifelse(Least_age>1,"adult", NA)))))

head(pres_stage)

pres_stage<-pres_stage %>% 
  mutate(date = as.Date(date))

#save this for later. Here we want to know if an individual is born in "other helgeland" or lurøy-onøy.
#this is so that we get the correct number of nestlings for each island.
born_locations<-pres_stage %>% 
  group_by(ID) %>% 
  arrange(Location, date) %>% 
  mutate(born_location = Location[which.min(date)]) %>% 
  mutate(born_island = Island[which.min(date)]) %>%
  ungroup() %>% 
  dplyr::select(ID, born_location, born_island) %>% 
  unique()



#subset years

pres_stage<-pres_stage %>% 
  filter(Year>=1994 & Year<=2022)

head(pres_stage)

#remove the islands we are not interested in.
pres_stage <-pres_stage %>% 
  filter(Island != 33) %>% 
  filter(Island != 0)


nrow(pres_clean)
nrow(pres_stage)

pres_stage %>% 
  dplyr::select(Location, Year) %>% 
  group_by(Location) %>% 
  summarise(min(Year), max(Year))


#In this analysis, individuals can change island, but they can only have one island
#per lifestage.
#First we make sure nestlings are only registered in one island
nest_floks<-pres_stage %>% 
  filter(lifestage == "nest") %>% 
  dplyr::select(ID,Year, Island, Location) %>% 
  unique()

nest_floks %>%
  group_by(ID) %>% 
  filter(n()>1) #five nestlings are registered in more than 1 island.
#remove double individuals
nest_floks<-nest_floks %>%
  group_by(ID) %>% 
  filter(n()==1) %>% 
  ungroup()
nrow(nest_floks)
length(unique(nest_floks$ID))

##adults
rec_ad_floks<-pres_stage %>% 
  filter(lifestage %in% c("recruit","adult")) %>% 
  dplyr::select(ID, Island, Location) %>% 
  unique()

length(unique(rec_ad_floks$ID))

rec_ad_floks %>%
  group_by(ID) %>% 
  filter(n()>1) #we have a few

rec_ad_floks<-pres_stage %>% 
  filter(lifestage %in% c("recruit","adult")) %>% 
  dplyr::select(ID,date,Year, Island, Location) %>% 
  unique()


nrow(rec_ad_floks)

double_island_same_year_adult<-rec_ad_floks %>%
  dplyr::select(-date) %>% 
  unique() %>% 
  group_by(ID, Year) %>% 
  filter(n()>1) %>% 
  ungroup()
length(unique(double_island_same_year_adult$ID))

rec_ad_floks<-rec_ad_floks %>% 
  dplyr::select(-date) %>% 
  unique()

##We will estimate true N inside Stan, so we want the observed number of individuals ####

rec_ad_floks %>% 
  group_by(Location) %>%
  summarise(min(Year))


rec_ad_floks<-rec_ad_floks %>% 
  group_by(Location, Year) %>% 
  mutate(N_obs = length(ID)) %>% 
  ungroup()



ages <- pres_stage %>% 
  dplyr::select(ID, Year, Least_age) %>% 
  unique()

rec_ad_floks<-rec_ad_floks %>% 
  left_join(ages)


#remove adults who are observed in more than two islands in a single year from the adult dataset
rec_ad_floks<-rec_ad_floks %>% 
  filter(!ID %in% double_island_same_year_adult$ID)
nrow(rec_ad_floks)
#remove the adult dispersers from the nest fledge datasets
fledge_floks<-fledge_floks %>% 
  filter(!ID %in% double_island_same_year_adult$ID)

nest_floks<-nest_floks %>% 
 filter(!ID %in% double_island_same_year_adult$ID) %>% 
  mutate(Least_age = 0)

#now we have removed double observations per year etc, we can left_join the born islands. 
#This is for calculating nestling production using surviving recruits
obs_rec<-rec_ad_floks %>% 
  left_join(born_locations)

#for observed recruits we get the island they were born in instead of where they were observed as recruits. 
#When we back-calculate nestlings we want the island where an individual was born, not where it dispersed.

obs_rec<-obs_rec %>% 
  group_by(Year, born_location) %>% 
  mutate(Rec_obs = length(ID[Least_age == 1])) %>% 
  ungroup() %>% 
  filter(!born_location %in% c("other helgeland", "lurøy-onøy")) 

obs_rec <- obs_rec %>% 
  dplyr::select(born_location, Year, Rec_obs) %>% 
  unique() %>% 
  rename(Location = born_location)%>% 
  arrange(Location, Year)

obs_pop<-rec_ad_floks %>% 
  dplyr::select(Location,Island, Year, N_obs) %>% 
  unique() %>% 
  arrange(Location, Year)


obs_pop<-obs_pop %>% 
  arrange(Location, Year)

obs_pop<-obs_pop %>% 
  left_join(obs_rec)



nest_floks$Least_age <- 0

head(rec_ad_floks)


rec_ad_floks %>% 
  dplyr::select(Location, Year) %>% 
  unique() %>% 
  group_by(Location) %>% 
  summarise(min(Year))

nest_rec_ad_floks<-rec_ad_floks %>% 
  bind_rows(nest_floks) %>% 
  arrange(Year)


nest_rec_ad_floks %>% 
  filter(is.na(N_obs))

stopifnot(
1>length(which(is.na(obs_pop$N_obs)))
)


double_test<-nest_rec_ad_floks %>% 
  group_by(Year, ID) %>% 
  filter(n()>1)
stopifnot(
nrow(double_test)==0
)



rec_ad_CH<-nest_rec_ad_floks %>% 
  dplyr::select(ID, Year) %>% 
  pivot_wider(names_from = Year, values_from = Year)

full_CH<-rec_ad_CH %>% 
  mutate(across(2:ncol(rec_ad_CH), ~as.numeric(ifelse(is.na(.x), "0", "1"))))

#nestling design matrix
nestling_year<-pres_stage %>% 
  dplyr::select(ID, Least_hatchyear) %>% 
  unique() %>% 
  filter(ID %in% nest_rec_ad_floks$ID) %>% 
  arrange(Least_hatchyear) %>%
  pivot_wider(values_from = Least_hatchyear, names_from = Least_hatchyear)

nestling_year<-nestling_year %>% 
  mutate(across(2:ncol(nestling_year), ~as.numeric(ifelse(is.na(.x), "0", "1"))))

#some individuals have hatchyears before the start of the CH. We cannot use the last year.
nestling_year<-nestling_year %>% 
  dplyr::select(ID,`1994`:`2021`)

##Adult age matrix ####

age<-nest_rec_ad_floks %>% 
  dplyr::select(ID,Least_age,Year) %>% 
  mutate(Least_hatchyear = Year - Least_age) %>% 
  pivot_wider(names_from = Year, values_from = Least_age) %>% 
  mutate(across(3:last_col(), ~ as.numeric(.)))

mean_adult_age<-mean(rec_ad_floks$Least_age)


age_filled<-age %>% 
  pivot_longer(cols = 3:last_col(), names_to = "Year", values_to = "age") %>%
  mutate(Year = as.numeric(Year)) %>% 
  group_by(ID) %>% 
  mutate(age = ifelse(
    Year>Least_hatchyear, Year - Least_hatchyear - mean_adult_age, #centre on the metapopulation mean adult age.
    NA)) %>% #individuals get NA when they are juveniles.
  filter(!is.na(Year)) %>%
  pivot_wider(names_from = Year, values_from = age) %>% 
  ungroup() %>% 
  mutate(across(3:last_col(), ~ replace_na(., 0)))

age_filled<-age_filled %>% 
  dplyr::select(ID,`1994`:`2021`)


##Adult design matrix####
adult_year<-pres_stage %>% 
  dplyr::select(ID, Least_hatchyear) %>% 
  unique() %>% 
  filter(ID %in% nest_rec_ad_floks$ID) %>% 
  arrange(Least_hatchyear) %>%
  mutate(Adult_year = Least_hatchyear +1) %>% 
  dplyr::select(ID, Adult_year) %>% 
  pivot_wider(values_from = Adult_year, names_from = Adult_year) %>% 
  pivot_longer(c(2:last_col()), names_to = "Year", values_to = "Adult_year") %>% 
  group_by(ID) %>% 
  fill(Adult_year, .direction = "down") %>%  
  ungroup() %>% 
  pivot_wider(values_from = Adult_year, names_from = Year)


adult_year<-adult_year %>% 
  mutate(across(2:ncol(adult_year), ~as.numeric(ifelse(is.na(.x), "0", "1"))))

adult_year<-adult_year %>% 
  dplyr::select(ID,`1994`:`2021`)

ncol(adult_year)
ncol(nestling_year)

colnames(pres_stage)

##Create sex column. Use genetically corrected sex for the years when we have it ####


sex_ecol<-pres_stage %>% 
  dplyr::select(ID, scriptsex) %>% 
  unique() %>% 
  mutate(sex2 = ifelse(scriptsex == "f", 0,
                      ifelse(scriptsex == "pf", 0,
                      ifelse(scriptsex == "m", 1,
                      ifelse(scriptsex == "pm", 1,NA)))))

sex_all <- sex_ecol %>% 
  left_join(sex_gen_corr)


sex_all<-sex_all %>% 
  group_by(ID) %>% 
  mutate(sex_combined = ifelse(!is.na(sex), sex, sex2)) %>% 
  dplyr::select(ID, sex_combined) %>% 
  ungroup() %>% 
  unique()

sex_all<-sex_all %>% 
  arrange(ID, sex_combined) %>% 
  group_by(ID) %>% 
  fill(sex_combined, .direction = "down") 
  
sex<-full_CH %>% 
  dplyr::select(ID) %>% 
  left_join(sex_all)

#check if all individuals who reach adulthood have a sex

max_ages<-pres_stage %>% 
  group_by(ID) %>% 
  mutate(max_age = max(Least_age)) %>% 
  ungroup() %>% 
  dplyr::select(ID, max_age) %>% 
  unique()

max_ages<- max_ages%>% 
  left_join(sex)

ads_no_sex<-max_ages %>% 
  filter(max_age>0) %>% 
  filter(ID %in% full_CH$ID) %>% 
  filter(is.na(sex_combined))

length(unique(ads_no_sex$ID)) #90 adults don't have sex. Assign a random sex

#90 adults don't have sex. Assign a random sex to those so that they can inform recapture rates
sex<-sex %>% 
  mutate(sex_combined = ifelse(ID %in% ads_no_sex$ID, rbinom(90,1,0.5), sex_combined))
sex %>% 
  filter(ID %in% ads_no_sex$ID) %>% 
  print(n=90)

#Since we cannot have NA in stan we need to assign a value to NA's in sex (all juveniles who never recruit). 
#Needs an adult design matrix to remove these in Stan.


sex$sex_combined[is.na(sex$sex_combined)] <- 0

##Island matrix. Individuals can change island
nest_rec_ad_floks$Island_stan <- as.numeric(as.factor(nest_rec_ad_floks$Island))

nest_rec_ad_floks %>% 
  dplyr::select(ID, Location) %>% 
  arrange(ID) %>%
  unique() %>% 
  group_by(ID) %>% 
  filter(n()>1) %>% 
  ungroup() %>% 
  dplyr::select(ID) %>% 
  unique() %>% 
  nrow()

length(unique(nest_rec_ad_floks$ID))

#we only have 534 dispersers.

islands_wide<-nest_rec_ad_floks %>% 
  dplyr::select(ID, Year, Island_stan) %>% 
  pivot_wider(names_from = Year, values_from = Island_stan)

#fill with the last known island, then fill with the first known island for individuals not caught as nestlings 
islands_filled<-islands_wide %>% 
  pivot_longer(2:ncol(islands_wide), names_to = "Year", values_to = "Island_stan") %>% 
  group_by(ID) %>% 
  fill(Island_stan, .direction="downup") %>% 
  ungroup()


islands_filled<-islands_filled %>%
  filter(Year<max(Year)) #remove the last year.we cannot assess survival here. But we can assess p. This means that the FY codes will mean different years for phi and p.

islands<-islands_filled %>% 
  pivot_wider(names_from = Year, values_from = Island_stan)

islands %>% 
  filter(`1994`!=`2021`) %>% 
  nrow() #without 2022 we only have 477 dispersers.

FY2<-islands_filled %>% 
  unite(FY, c(Island_stan, Year), remove = F) %>% 
  mutate(FY = as.numeric(as.factor(FY))) %>% 
  dplyr::select(ID,Year, FY) %>%
 pivot_wider(names_from = Year, values_from = FY)

##density matrix

obs_pop$Island_stan <- as.numeric(as.factor(obs_pop$Island))

demo4<-obs_pop %>% 
  dplyr::select(Island, Location,Island_stan,Year, N_obs, Rec_obs) %>% 
  unique()


demo4<-demo4 %>% 
  unite("Island_stan_Year", c("Island_stan","Year"), remove = F)

demo4$FY <- as.numeric(as.factor(demo4$Island_stan_Year))


N<-demo4 %>% 
  dplyr::select(Island_stan, N_obs, Year) %>% 
  unique() %>% 
  arrange(Year) %>% 
  filter(Year<max(Year)) %>% #remove the last year. We do not use that column since we cannot assess survival then.
  pivot_wider(names_from = Year, values_from = N_obs) %>% 
  arrange(Island_stan)
#change NA to 1. Most of these values will be excluded in the loop.
N[is.na(N)] <- 1
N[N == 0] <- 1 #change 0 to 1. Most of these values will be excluded in the loop.


N_rec<-demo4 %>% 
  dplyr::select(Island_stan, Rec_obs, Year) %>% 
  unique() %>% 
  arrange(Year) %>% 
  filter(Year<max(Year)) %>% #remove the last year. We do not use that column since we cannot assess survival then.
  pivot_wider(names_from = Year, values_from = Rec_obs) %>% 
  arrange(Island_stan)
#change NA to 1. Most of these values will be excluded in the loop.
N_rec[is.na(N_rec)] <- 1
N_rec[N_rec == 0] <- 1 #change 0 to 1. Most of these values will be excluded in the loop.

unique(demo4$Location)
length(unique(demo4$Location))

inout<-demo4 %>% 
  mutate(in_out = ifelse(Location %in% c("aldra", "gjerøy", "hestmannøy", "indre kvarøy", "nesøy"),1,0)) %>% 
  dplyr::select(Island,Island_stan, in_out) %>% 
  unique() %>% 
  arrange(Island_stan)

#Start year for each island (when nestling surveys started)

start_year<-nest_rec_ad_floks %>% 
  dplyr::select(Location, Island_stan, Year) %>% 
  unique() %>% 
  mutate(Year_stan = as.numeric(as.factor(Year))) %>% 
  group_by(Location,Island_stan) %>% 
  summarise(start_year_stan = min(Year_stan)) %>% 
  ungroup() %>% 
  arrange(Island_stan)
  
##create stan data ####
stopifnot(
identical(full_CH$ID, islands$ID)
)


nest_mat <- nestling_year[match(full_CH$ID, nestling_year$ID), ]
adult_mat <- adult_year[match(full_CH$ID, adult_year$ID), ]
age_mat <- age_filled[match(full_CH$ID, age_filled$ID), ]
sex <- sex[match(full_CH$ID, sex$ID), ]


FY2 <- FY2[match(full_CH$ID, FY2$ID), ]

stopifnot(c(
identical(nest_mat$ID, full_CH$ID),
identical(age_mat$ID, full_CH$ID),
identical(sex$ID, full_CH$ID),
identical(adult_mat$ID, full_CH$ID),
identical(FY2$ID, full_CH$ID),
identical(N$Island_stan, inout$Island_stan),
identical(N$Island_stan, start_year$Island_stan)
)
)

CH<-full_CH %>% 
  dplyr::select(-ID)

FY2<-FY2 %>% 
  dplyr::select(-ID)
nFY2<-max(FY2)

islands<-islands %>% 
  dplyr::select(-ID)

island_id <-1:max(islands)

islands2<-as.data.frame(islands)[,1]

start_year<-start_year %>% 
  dplyr::select(-Location)

nest_mat<-nest_mat %>% 
  dplyr::select(-ID)

adult_mat<-adult_mat %>% 
  dplyr::select(-ID)

age_mat<-age_mat %>% 
  dplyr::select(-ID)

age_sq_mat <- age_mat^2

N_stan<-N %>% 
  dplyr::select(-Island_stan)

N_rec<-N_rec %>% 
  dplyr::select(-Island_stan)
nFY<-nrow(N_stan)*ncol(N_stan)
FY <- as.data.frame(N_stan)
FY[]<-1:nFY
colnames(FY)<-NULL #Just to avoid confusion.Naming doesn't matter. FY references different years for FY_phi and FY_p. Recapture is assessed up to the last year but survival is assessed up to the last year - 1.

nFY<-nrow(FY)*ncol(FY) #remember! The codes won't point to the same year for p and phi. phi excludes the last year and p the first.

stan_data <- list(y =CH,  nind = dim(CH)[1], n_occasions =
                    dim(CH)[2], flok=islands, nflok=max(islands), start_year = start_year$start_year_stan,
                  max_startyear = max(start_year$start_year_stan),
                  obs_pop=N_stan,obs_rec = N_rec,in_out=inout$in_out, 
                  sex = sex$sex_combined, age = age_mat, age_sq = age_sq_mat,
                  NR=nest_mat, AD=adult_mat, FY=FY, nFY=nFY 
)


##Run Stan model####


inits <- function() list(
  sigmaF_phi =  array(runif(2, 0, 1), dim = 2),
  sigmaF_phi_rec =  array(runif(2, 0, 1), dim = 2),
  sigmaFY_phi =  array(runif(1, 0, 1), dim = 1),
  sigmaFY_phi_rec =  array(runif(1, 0, 1), dim = 1),
  sigmaFY_p =  array(runif(1, 0, 1), dim = 1)) #  need to specify an explicit dimension using array and dim because we specify sigmas as vectors of dimensions 1 or 2! We can change sigmas to real and then we don't have to do it like this.

init_list <- inits()
str(init_list)

mod<-stan_model("Stan scripts/CMR nest to ad dispersers ELASTICITY gen quant block.stan")

params <- c("mu_phi","mu_phi_rec","B_sex","B_age", "B_age_sq","B_sex_age", "B_in_out", "B_in_out_rec", "gamma","gamma_rec","gamma_in_out", "gamma_in_out_rec", "mu_p", "Sigma2_YI_p",
            "Sigma2_isl_phi","Sigma2_isl_phi_rec","Sigma2_isl_gamma", "Sigma2_isl_gamma_rec", 
            "Sigma2_FY_phi", "Sigma2_FY_phi_rec", "cov_isl", "cov_isl_rec", "f", "f_simple", "s_n", "s_a",
            "elasticity","lambda","Nest_corr", "Nest_corr_F", "Rec_corr",
            "epsilonF_phi", "epsilonF_phi_rec","epsilonFY_phi", "epsilonFY_phi_rec","LF","LF_rec",  "zNR", "N_corr", "p")

ni <- 4000
nt <- 2
nc <- 4


stopifnot(
  nrow(CH)>15000
)

nest_to_rec_IPM_genquant_block<- sampling(mod,data = stan_data, pars = params,
                                          chains = nc, iter = ni, thin = nt, init = inits,
                                          cores=nc)
##OBS! Island Stan is created using the island column, not the location column like in the fledge script#####

saveRDS(nest_to_rec_IPM_genquant_block, "Workspace backup/nest_to_rec_IPM_genquant_block.rds")
nest_to_rec_IPM_genquant_block<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")
round(summary(nest_to_rec_IPM_genquant_block)$summary[c(1:45), c(6,4,8,9,10)],5)
rstan::traceplot(nest_to_rec_IPM_genquant_block, pars = c("cov_isl", "cov_isl_rec"))
launch_shinystan(nest_to_rec_IPM_genquant_block)
rstan::check_hmc_diagnostics(nest_to_rec_IPM_genquant_block)

library(bayesplot)
library(mcmcplots)
mcmc<-As.mcmc.list(nest_to_rec_IPM_genquant_block)
rstan::traceplot(nest_to_rec_IPM_genquant_block, pars =params[1:13])
denplot(mcmc, parms = c("Sigma2_isl_gamma_rec", "Sigma2_isl_phi_rec", "cov_isl_rec"))

nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Sigma2_isl_phi_rec, Sigma2_isl_gamma_rec, cov_isl_rec) %>%
  mutate(corr_rec = cov_isl_rec/(sqrt(Sigma2_isl_phi_rec)*sqrt(Sigma2_isl_gamma_rec))) %>% 
  median_qi(corr_rec)

nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Sigma2_isl_phi, Sigma2_isl_gamma, cov_isl) %>%
  mutate(corr_ad = cov_isl/(sqrt(Sigma2_isl_phi)*sqrt(Sigma2_isl_gamma))) %>% 
  median_qi(corr_ad)

x<-as.data.frame(round(summary(nest_to_rec_IPM_genquant_block)$summary[, c(6,4,8,9,10)],5))

x %>% 
  filter(is.na(Rhat))

x %>% 
  arrange(n_eff)
x %>% 
  arrange(Rhat)

##compare elasticities from stan with manually computed ####
elast_frame<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(elasticity[Island_stan, row, col]) %>% 
  group_by(Island_stan ,row, col) %>% 
  median_qi(elasticity) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, row, col, elasticity)

isl_loc<-demo4 %>% 
  dplyr::select(Island, Location) %>% 
  unique()

inout<-inout %>% 
  left_join(isl_loc)

elast_frame<-elast_frame %>% 
  left_join(inout)

elast_frame_sum<-elast_frame %>% 
  pivot_wider(names_from = c(row, col), values_from = elasticity) %>% 
  mutate(juv_sum = `1_1` + `1_2` +`2_1`) %>% 
  mutate(elasticity_ratio = juv_sum/`2_2`)

elast_frame<-elast_frame %>% 
  filter((row == 1 & col == 2) | (row == 2 & col == 2)) %>% 
  dplyr::select(-col) %>% 
  pivot_wider(names_from = row, values_from = elasticity) %>% 
  rename(Fec_juv = `1`) %>% 
  rename(Ad = `2`) %>% 
  mutate(elasticity_ratio = Fec_juv/Ad)

elast_frame
elast_frame_sum

f_simple<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(f_simple[Island_stan]) %>% 
  dplyr::select(Island_stan, f_simple) %>% 
  ungroup()
class(f_simple)
f_simple %>%
  median_qi(f_simple, .width =0.9)

f_simple<-f_simple %>% 
  left_join(inout)

f_simple %>% 
  group_by(in_out) %>% 
  median_qi(f_simple, .width =0.9)

f_simple<-f_simple %>% 
  group_by(Island_stan) %>% 
  median_qi(f_simple, .width =0.9) %>% 
  ungroup() %>% 
  arrange(f_simple)
f_simple


s_n<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(s_n[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(s_n) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, s_n)

s_a<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(s_a[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(s_a) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, s_a)

f_simple<-f_simple %>% 
  left_join(s_n) %>% 
  left_join(s_a)
f_simple<-f_simple %>% 
  left_join(inout)
f_simple %>% 
  arrange(in_out)
pm_list<-list()
for(i in 1:10){
  proj_mtrx<-matrix(c(0,  f_simple$f_simple[i], f_simple$s_n[i], f_simple$s_a[i]), nrow=2, byrow=TRUE,
                    dimnames=list(c("nestling","adult"),
                                  c("nestling","adult")))
  pm_list[[i]] <- proj_mtrx
  names(pm_list)[[i]] <- f_simple$Location[i]
}
library(popbio)

elasticity_list<-lapply(pm_list, elasticity)
elasticity_list
elast_frame


#check blups
nest_to_rec_IPM_genquant_block %>% 
  spread_draws(epsilonF_phi_rec[variable, Island_stan]) %>% 
  filter(variable == 2) %>% 
  group_by(Island_stan) %>% 
  median_qi(epsilonF_phi_rec) %>% 
  left_join(inout)


nest_to_rec_IPM_genquant_block %>% 
  spread_draws(epsilonF_phi[variable, Island_stan]) %>% 
  filter(variable == 2) %>% 
  group_by(Island_stan) %>% 
  median_qi(epsilonF_phi) %>% 
  left_join(inout)

gammas<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(gamma,gamma_rec,gamma_in_out, gamma_in_out_rec, epsilonF_phi[variable, Island_stan], epsilonF_phi_rec[variable, Island_stan]) %>% 
  filter(variable == 2)

gammas<-gammas %>% 
  left_join(inout)

#with an adult design matrix gamma_rec and epsilonF_phi_rec are not deviations
gamma_ratios<-gammas %>% 
  mutate(gamma_ad_isl = ifelse(in_out == 0, gamma + epsilonF_phi, gamma  +gamma_in_out+ epsilonF_phi)) %>% 
  mutate(gamma_nest_to_rec_isl = ifelse(in_out == 0,gamma_rec + epsilonF_phi_rec, gamma_rec +gamma_in_out_rec + epsilonF_phi_rec)) %>% 
  mutate(gamma_nest_ad_ratio = gamma_nest_to_rec_isl/gamma_ad_isl) %>% 
  group_by(Island_stan) %>% 
  median_qi(gamma_nest_ad_ratio,gamma_ad_isl,gamma_nest_to_rec_isl, .width = 0.9) %>% 
  dplyr::select(Island_stan, gamma_nest_ad_ratio, gamma_ad_isl, gamma_nest_to_rec_isl)

gamma_ratios<-gamma_ratios %>% 
  left_join(inout)

gamma_ratios

elast_frame

elast_frame_sum<-elast_frame_sum %>% 
  dplyr::select(Location, elasticity_ratio) %>% 
  rename(elasticity_ratio_sum = elasticity_ratio)

elast_ratio <- elast_frame %>% 
  left_join(gamma_ratios)

elast_ratio<-elast_ratio %>% 
  left_join(elast_frame_sum)

elast_ratio<-elast_ratio %>% 
  arrange(in_out)

rec_isl_gamma<-readRDS("Workspace backup/rec_isl_gamma.rds")

rec_isl_gamma<-rec_isl_gamma %>% 
  dplyr::select(Location, gamma_fledge_to_rec_isl)

elast_ratio <- elast_ratio %>% 
  left_join(rec_isl_gamma)

elast_ratio<-elast_ratio %>% 
  mutate(gamma_fledge_to_rec_ad_ratio = gamma_fledge_to_rec_isl/gamma_ad_isl)

##Elasticity figures ####
par(mfrow =c(1,1))
par(mar=c(4,4,4,4))

pdf(file = "Figures empirical data/Nest to rec adult gamma elasticity ratio new N.pdf",   # The directory you want to save the file in
    width = 6.5, # The width of the plot in inches 6.43
    height = 5) # 3.92

plot(gamma_nest_ad_ratio  ~ log(elasticity_ratio),
     ylab ="",
     xlab="",
     data = elast_ratio, col = in_out+1, pch=16)

mtext(expression(paste("Strength of density dep. ratio ",
                       frac(gamma[nr], gamma[ad]))), 
      side = 2, line = 1.7, cex = 1)

mtext(expression(paste("Elasticity ratio Log", bgroup("(",frac(e('F'), e(phi[ad])), ")"))), #we can call this F because fecundity and juvenile survival have the same elasticities. F is therefore recruit production including fecundity.
      side = 1, line = 4, cex = 1.1)  # Increase line number to move it further down
library(basicPlotteR)
addTextLabels(log(elast_ratio$elasticity_ratio) , elast_ratio$gamma_nest_ad_ratio , str_to_title(elast_ratio$Location),col.label = elast_ratio$in_out+1,
              col.line = "black")
legend("topright", legend=c("Inner farm islands", "Outer non-farm islands"),
       col=c("red", "black"), pch = 16, cex=1)
dev.off()

pdf(file = "Figures empirical data/Fledge to rec adult gamma elasticity ratio new N.pdf",   # The directory you want to save the file in
    width = 6.5, # The width of the plot in inches 6.43
    height = 5) # 3.92

plot(log(gamma_fledge_to_rec_ad_ratio)  ~ elasticity_ratio, data = elast_ratio,
     ylab = "", 
     xlab = "",
     col = in_out+1, pch=16)

mtext(expression(paste("Strength of density dep. ratio Log",
                       bgroup("(",frac(gamma[fr], gamma[ad]), ")"))), 
      side = 2, line = 1.7, cex = 1)

mtext(expression(paste("Elasticity ratio Log", bgroup("(",frac(e('F'), e(phi[ad])), ")"))), #we can call this F because fecundity and juvenile survival have the same elasticities. F is therefore recruit production including fecundity.
      side = 1, line = 4, cex = 1.1)  # Increase line number to move it further down

addTextLabels(log(elast_ratio$elasticity_ratio) , log(elast_ratio$gamma_fledge_to_rec_ad_ratio) , str_to_title(elast_ratio$Location),col.label = elast_ratio$in_out+1,
              col.line = "black")
legend("topright", legend=c("Inner farm islands", "Outer non-farm islands"),
       col=c("red", "black"), pch = 16, cex=1)

dev.off()
#N_corr

nest_to_rec_IPM_genquant_block %>% 
  spread_draws(p[Island_stan, Year_stan]) %>% 
  group_by(Island_stan, Year_stan) %>% 
  median_qi(p)

N_corr<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(N_corr[Island_stan, Year_stan]) %>% 
  group_by(Island_stan, Year_stan) %>% 
  median_qi(N_corr, .width = 0.9) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, Year_stan, N_corr, .lower, .upper, .width)

N_corr


year_island_IDs<-N %>% 
  pivot_longer(c(2:last_col()), names_to = "Year", values_to = "N_obs") %>% 
  mutate(Year = as.numeric(Year)) %>% 
  mutate(Year_stan = as.numeric(as.factor(Year)))

isl_ids<-data_observed %>% 
  dplyr::select(Location, Island, in_out) %>% 
  filter(Island != 33) %>% 
  unique() %>% 
  mutate(Island_stan = as.numeric(as.factor(Island)))

year_island_IDs <- year_island_IDs %>% 
  left_join(isl_ids)


N_corr<-N_corr %>% 
  left_join(year_island_IDs)
demo_name<-demo %>% 
  rename(Location = Islandname)

demo_name %>% 
  filter(Location == "myken" & Year == 2014)

N_corr<-N_corr %>% 
  left_join(demo_name)

N_corr %>% 
  filter(Island_stan == 2) %>% 
  print(n =28)

N_corr_print <- N_corr %>% 
  dplyr::select(Location,Island, Year, N_corr, .lower, .upper, .width, N_obs)

N_corr_print[is.na(N_corr_print$Location),]

N_corr_print %>% 
  filter(Location == "aldra") %>% 
  print(n = 28)

N_corr_print %>% 
  filter(is.na(Location))
##Print estimated pop sizes ####

write_delim(N_corr_print, "Cleaned data for density dependence analyses/estimated_pop_sizes_nestling_to_adult_model.txt", delim = ";")
