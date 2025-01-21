#Adult CMR

library(rstan)
library(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
cmr_raw<-read.table("Data/cleaned_density_dependence_data_CMR_1998_2020.txt", sep=";", header=TRUE)


#Remove second year dispersers.
double_islands<-cmr_raw %>% 
  dplyr::select(ID, Location) %>% 
  unique() %>% 
  group_by(ID) %>% 
  filter(n()>1)

#also remove the filled in years (no field observations these years).
cmr_sub<-cmr_raw%>%
  filter(!(ID%in% double_islands$ID)) %>% 
  filter(filled=="no") %>% 
  arrange(Location,Year) %>% 
  unique()


age<-cmr_sub %>% 
  dplyr::select(ID, Least_hatchyear,Location,sex,Least_age,Year) %>% 
  pivot_wider(names_from = Year, values_from = Least_age) %>% 
  mutate(across(5:last_col(), ~ as.numeric(.)))

#Convert the ages and the NA's to 1's and 0's to generate a capture history. 
# Remove ID's (a row is an individual) and Location.
CH<-age %>% 
  mutate(across(5:last_col(), ~ replace_na(.,0))) %>% 
  mutate(across(5:last_col(), ~ ifelse(.x>0, 1, 0))) %>% 
  dplyr::select(-(1:2))



#Fill the age matrix from the first occasion to the last. 
# These are possible ages that an individual could have 
# given the first time it was caught. It is not realised age.

age_filled<-age %>% 
  pivot_longer(cols = 5:last_col(), names_to = "Year", values_to = "age") %>%
  mutate(Year = as.numeric(Year)) %>% 
  group_by(ID) %>% 
  mutate(age = ifelse(
    Year>Least_hatchyear, Year - Least_hatchyear-1, # so that the recruit year is 0. Then the intercept will be recruit survival
    NA)) %>% 
  filter(!is.na(Year)) %>%
  pivot_wider(names_from = Year, values_from = age) %>% 
  ungroup() %>% 
  mutate(across(5:last_col(), ~ replace_na(., 0))) %>% 
  dplyr::select(-last_col()) # remove the last year. We do not use that column since we cannot assess survival for the last year


max(age_filled$`1998`) #we have a bunch of individuals who enter the data when they are old.

#convert to matrix
age_matrix<-age_filled %>% 
  dplyr::select(5:last_col()) %>% 
  as.matrix()

#square age
age_matrix_sq <- age_matrix^2 

#Create a pop size matrix with islands on the rows and years as cols
#We need the density at t, because the survival is assessed for the 
#transition after an individual was captured.


#Also, mean-centred density needs to be recalculated for the years used.
#We only use 2020 for recapture, not survival.

popsize<-cmr_raw %>% 
  dplyr::select(Location, Year, N) %>%
  filter(Year<max(cmr_raw$Year)) %>% # remove the last year. We do not use that column since we cannot assess survival for the last year
  unique() %>% 
  arrange(Location, Year) #sort the dataframe so that islands (the rows) are in ascending order (starting with Aldra), Years are also in ascending order

popsize<-popsize %>% 
  group_by(Location) %>% 
  mutate(Mean_N_isl = mean(N)) %>% 
  ungroup() %>% 
  mutate(mc_N = N - Mean_N_isl)
mean(popsize$mc_N) #basically 0

#scale density to 1 sd (rounded to 40)
scaler = 40
zNR_long<-popsize %>% 
  mutate(zNR = mc_N/scaler) 


#Create popsize matrix
pop_year_island<-zNR_long %>% 
  dplyr::select(Location, Year, zNR) %>% 
  pivot_wider(names_from = Year, values_from = zNR) %>% 
  dplyr::select(-Location) %>% 
  as.matrix()

head(CH)
nrow(pop_year_island)
nrow(age_matrix)

#create island and sex vectors
island<-as.numeric(as.factor(CH$Location)) #save Locations as vector
sex <-as.numeric(CH$sex)

#Create inner_outer vector.
ios<-cmr_sub %>% 
  dplyr::select(Location, in_out) %>% 
  unique()
io <- CH %>% 
  left_join(ios, by = "Location") %>% 
  dplyr::select(in_out)


CH <- CH[,3:ncol(CH)]#and now remove the Location and sex columns

#Create a matrix with unique numbers for each island-year combination
nFY<-nrow(pop_year_island)*ncol(pop_year_island)

FY<-pop_year_island
FY[]<-1:nFY
colnames(FY)<-NULL #Just to avoid confusion.Naming doesn't matter. FY is made from the pop size frame. It references different years for FY_phi and FY_p. Recapture is assessed up to the last year (2020) but survival is assessed up to 2019.

FY
sex
island
nFY
dim(CH)[2]
colnames(CH)
dim(FY)[2]
nrow(CH)
dim(pop_year_island)
head(pop_year_island)
head(CH)
head(age_matrix)
head(age_matrix_sq)
nrow(age)
nrow(io)


##Full analysis ####



stan_data <- list(y =CH,  nind = dim(CH)[1], n_occasions =
                    dim(CH)[2], flok=island, nflok=max(island), age=age_matrix,age_sq=age_matrix_sq, sex=sex, 
                  n=pop_year_island,in_out=io$in_out, FY=FY, nFY=nFY 
)

#Parameters monitored

params <- c("mu_phi","mu_p", "gamma","B_age", "B_age_sq", "B_sex", "B_in_out","interac_in_out_gamma", "interac_age_sex","Sigma2_Isl_phi","Sigma2_Isl_gamma", "Sigma2_YI_phi","Sigma2_YI_p", "cov_Isl", "isl_phi", "isl_g","p")

## MCMC settings

ni <- 12000
nt <- 30
nc <- 4

mod<-stan_model(file ="Stan scripts/CMR density reg adult surv Paul_old_no_ind_int.stan")


## Call Stan from R 
cjs_temp_raneff2 <- sampling(mod,
                             data = stan_data, pars = params,
                             chains = nc, iter = ni, thin = nt,
                             cores=nc)

saveRDS(cjs_temp_raneff2, file("Workspace backup/CMR_adults_interacs_old_Paul_age_sq_zNR.rds")) #The "old" version is the correct one.
