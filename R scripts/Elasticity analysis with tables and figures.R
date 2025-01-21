#Elasticities

library(tidyverse)
library(rstan)
library(tidybayes)
library(popbio)
library(arm)
library(readxl)

source("Functions/CMR2table_juveniles.R")
paramlatex<-read_xlsx(path="Tables empirical data/Parameter to latex.xlsx")

CMR_adult_survival<-readRDS("Workspace backup/CMR_adults_interacs_old_Paul_age_sq_zNR.rds")
CMR_juv<-readRDS("Workspace backup/CMR_juvenile_survival_Tor_zNR.rds")
data_observed<-read.table("Data/cleaned_density_dependence_data.txt", sep=";", header=TRUE)
#get island blups

#for juveniles I sorted the islands by the numerical Island column for the stan data.
isl_sys<-data_observed %>% 
  filter(Location != "lurøy-onøy") %>% 
  dplyr::select(Island, Location, in_out) %>% 
  arrange(Island) %>% 
  unique()

isl_sys$islcode <- 1:nrow(isl_sys)
CMR_juv %>% 
  get_variables()

#I forgot to add inner outer to the BLUPs in stan so I do that here.
juv_draws<-CMR_juv %>% 
  spread_draws(phi_fledge_interac_inout,isl_phi_fledge[islcode],phi_rec_interac_inout, isl_phi_rec[islcode])
juv_draws<-juv_draws %>% 
  left_join(isl_sys, by = "islcode")

nestling_surv<-juv_draws %>% 
  group_by(Location) %>% 
  mutate(isl_blups = ifelse(in_out == 0, isl_phi_fledge,
                            phi_fledge_interac_inout + isl_phi_fledge)) %>% 
  median_qi(isl_surv = invlogit(isl_blups), .width = 0.9) %>% 
  as.data.frame()

fledgling_surv<-juv_draws %>% 
  group_by(Location) %>% 
  mutate(isl_blups = ifelse(in_out == 0, isl_phi_rec,
                            phi_rec_interac_inout + isl_phi_rec)) %>% 
  median_qi(isl_surv = invlogit(isl_blups), .width = 0.9) %>% 
  as.data.frame()

#For adults I sorted the islands by the character column Location.
chr_locs<-data_observed %>% 
  dplyr::select(Location, in_out) %>% 
  arrange(Location) %>% 
  unique()

chr_locs$islcode <- 1:nrow(chr_locs)

CMR_adult_survival %>% 
  get_variables()
#I forgot to add inner_outer to the BLUPS in stan
ad_draws<-CMR_adult_survival %>% 
  spread_draws(isl_phi[islcode], B_in_out,  B_age, B_age_sq, B_sex, interac_age_sex)

ad_draws<-ad_draws %>% 
  left_join(chr_locs, by = "islcode")

cmr_raw<-read.table("Data/cleaned_density_dependence_data_CMR_1998_2020.txt", sep=";", header=TRUE)
mean_ages<-cmr_raw %>% 
  group_by(in_out) %>% 
  summarise(age = mean(Least_age)) %>% 
  ungroup()

inner_age<-mean_ages[2,2]$age
outer_age<-mean_ages[1,2]$age

ad_surv<-ad_draws %>% 
  group_by(Location) %>% 
  mutate(isl_blups = ifelse(in_out == 0, isl_phi + B_age*outer_age + B_age_sq*outer_age^2 + B_sex*0.5 + interac_age_sex*0.5*outer_age,
                            B_in_out + isl_phi + B_age*inner_age + B_age_sq*inner_age^2 + B_sex*0.5 + interac_age_sex*0.5*inner_age)) %>% 
  median_qi(isl_surv = invlogit(isl_blups), .width = 0.9)

#remove lurøy onøy since we don't have juveniles there
ad_surv<-ad_surv %>% 
  filter(Location != "lurøy-onøy") %>% 
  as.data.frame()

###Back-calculate number of nestlings####
stage_location_year_summary<-read.table("Data/stage_location_year_summary.txt",sep=";", header = T, stringsAsFactors = F)
growth_factors <-read.table("Data/growth_factors.txt",sep=";", header = T, stringsAsFactors = F)
unique(growth_factors$Location)

FY<-stage_location_year_summary %>% 
  dplyr::select(Location, Year) %>% 
  unique() %>% 
  arrange(Location, Year) %>% 
  unite(FY, c(Location, Year))

FY2<-growth_factors %>% 
  dplyr::select(Location, Year) %>% 
  unique() %>% 
  arrange(Location, Year)%>% 
  unite(FY, c(Location, Year))

stopifnot(
identical(FY$FY, FY2$FY)
)

stage_location_summary<-stage_location_year_summary %>% 
  group_by(stage, Location) %>% 
  summarise(total =sum(total)) %>% 
  ungroup()

reverse_inner_outer = "no"
add_difference_to_intercept = "yes"
width = 0.9
juv_res<-CMR2table_juveniles(CMR_juv,paramlatex,reverse_inner_outer,add_difference_to_intercept, width)

(nest2fledge_res<-juv_res[[1]])


#here we need the island BLUPs for recapture probability
CMR_juv_p<-readRDS("Workspace backup/CMR_juvenile_survival_Tor_zNR_tracked_p.rds")
FY_identifiers_juvenile_CMR<-readRDS("Workspace backup/FY_identifiers_juvenile_CMR.rds")

FY_identifiers_juvenile_CMR<-FY_identifiers_juvenile_CMR %>%
  separate(FY, c("Island", "Year"))

isl_loc<-data_observed %>% 
  dplyr::select(Island, Location) %>%
  mutate(Island = as.character(Island)) %>% 
  unique()

FY_identifiers_juvenile_CMR<-FY_identifiers_juvenile_CMR %>% 
  left_join(isl_loc) %>% 
  arrange(Island, Year)

fledge_p<-CMR_juv_p %>% 
  spread_draws(mu_p, epsilonFY_p_fledge[FY_code])
fledge_p<-fledge_p %>% 
  left_join(FY_identifiers_juvenile_CMR, by = "FY_code")


isl_fledge_p<-fledge_p %>% 
  group_by(Location) %>% 
  mutate(isl_fledge_p = mu_p + epsilonFY_p_fledge) %>% 
  median_qi(isl_fledge_p = invlogit(isl_fledge_p), .width = 0.9) %>% 
  ungroup() %>% 
  dplyr::select(Location, isl_fledge_p)


location_fledge <-stage_location_summary %>% 
  filter(stage == "fledge")

location_fledge<-location_fledge %>% 
  left_join(isl_fledge_p)

location_fledge<-location_fledge %>% 
  mutate(real_fledglings = total/isl_fledge_p)

nest_to_fledge_loc<-nestling_surv %>% 
  dplyr::select(Location, isl_surv)

location_fledge_nest<-location_fledge %>% 
  left_join(nest_to_fledge_loc) %>% 
  mutate(real_nestlings = real_fledglings/isl_surv)

adults_tot_loc<-growth_factors %>% 
  group_by(Location) %>% 
  summarise(tot_adults = sum(N)) %>% 
  ungroup()%>% 
  as.data.frame()

location_fledge_nest<-location_fledge_nest %>% 
  left_join(adults_tot_loc) %>% 
  mutate(fecundity = real_nestlings/tot_adults)%>% 
  as.data.frame()

fecundity_locations<- location_fledge_nest %>% 
  dplyr::select(Location, fecundity) %>% 
  as.data.frame()
str(fecundity_locations)
##Create projection matrices####
#loop over islands
pm_list<-list()
fecundity_locations[1,2]
nestling_surv
fledgling_surv
ad_surv[1,2]
for(i in 1:10){
  proj_mtrx<-matrix(c(0, 0, fecundity_locations[i,2], nestling_surv[i,2], 0, 0, 0, fledgling_surv[i,2], ad_surv[i,2]), nrow=3, byrow=TRUE,
                    dimnames=list(c("nestling","fledgling","adult"),
                                  c("nestling","fledgling","adult")))
  pm_list[[i]] <- proj_mtrx
  names(pm_list)[[i]] <- fecundity_locations[i,1]
}
library(popbio)

elasticity_list<-lapply(pm_list, elasticity)
sensitivity_list<-lapply(pm_list, sensitivity)
lambda_frame<- data.frame(lambda= sapply(pm_list, lambda))
lambda_frame

#extract adult island specific gammas
ad_gamma_draws<-CMR_adult_survival %>% 
  spread_draws(interac_in_out_gamma, isl_g[islcode])

ad_gamma_draws<-ad_gamma_draws%>% 
  left_join(chr_locs, by = "islcode")

ad_gammas<-ad_gamma_draws %>% 
  group_by(Location) %>% 
  mutate(gamma_blup =ifelse(in_out == 0, isl_g,
                            interac_in_out_gamma + isl_g)) %>% 
  median_qi(island_gamma = gamma_blup, .width =0.9) %>% 
  ungroup() %>% 
  filter(Location != "lurøy-onøy")
elast_ratio<-matrix(NA, 10, 4)
colnames(elast_ratio) <- c("Location", "F_Ad_ratio", "Rec_Ad_ratio", "Fledge_ad_ratio")
elast_ratio[,1]<-ad_gammas$Location

for(i in 1:10){
  elast_ratio[i,2]<-elasticity_list[[i]][1,3]/elasticity_list[[i]][3,3]
  elast_ratio[i,3]<-elasticity_list[[i]][3,2]/elasticity_list[[i]][3,3]
  elast_ratio[i,4]<-elasticity_list[[i]][2,1]/elasticity_list[[i]][3,3]
}

elast_ratio<-as.data.frame(elast_ratio)
elast_ratio[,2]<-as.numeric(elast_ratio[,2])
elast_ratio[,3]<-as.numeric(elast_ratio[,3])
elast_ratio[,4]<-as.numeric(elast_ratio[,4])

ad_gammas<-ad_gammas %>% 
  left_join(elast_ratio) %>% 
  left_join(chr_locs) %>% 
  dplyr::select(-.lower, -.upper, -.width, -.point, -.interval)

#Extract fledgling to recruit gammas

rec_gamma_draws<-CMR_juv %>% 
  spread_draws(isl_g_rec[islcode], gamma_rec_interac_inout)

rec_gamma_draws<-rec_gamma_draws %>% 
  left_join(isl_sys, by = "islcode")

rec_gammas<-rec_gamma_draws %>% 
  group_by(Location) %>% 
  mutate(isl_rec_gamma = ifelse(in_out == 0, isl_g_rec,
                                gamma_rec_interac_inout  +isl_g_rec)) %>% 
  median_qi(isl_rec_gamma, .width=0.9) %>% 
  ungroup() %>% 
  dplyr::select(Location, isl_rec_gamma)

ad_rec_gammas <- ad_gammas %>% 
  left_join(rec_gammas)

ad_rec_gammas<-ad_rec_gammas %>% 
  mutate(gamma_ratio =  isl_rec_gamma/island_gamma)

##Figure elasticity ratio ####
#We expect a strong (more negative) gamma with a high ratio.
#F and juvenile survival have the same elasticities
par(mfrow = c(1,1))
par(mar= c(5,5,4,4))
pdf(file = "Figures empirical data/Gamma elasticity ratio.pdf",   # The directory you want to save the file in
    width = 6.43, # The width of the plot in inches
    height = 3.92)
plot(log(gamma_ratio)~log(Rec_Ad_ratio), data = ad_rec_gammas, 
     ylab = "", 
     xlab = "",
     cex.lab = 1,
     col = in_out+1,
     pch = 16)

mtext(expression(paste("Strength of density dep. ratio Log",
                       bgroup("(",frac(gamma[fl], gamma[ad]), ")"))), 
      side = 2, line = 1.7, cex = 1)

mtext(expression(paste("Elasticity ratio Log", bgroup("(",frac(e('F'), e(P[ad])), ")"))), #we can call this F because fecundity and juvenile survival have the same elasticities. F is therefore recruit production including fecundity.
      side = 1, line = 4, cex = 1.1)  # Increase line number to move it further down

legend(0.15, 0.7, legend=c("Inner islands", "Outer islands"),
       col=c("red", "black"), pch = 16, cex=1)

abline(lm(log(gamma_ratio)~log(Rec_Ad_ratio), data = ad_rec_gammas), lwd = 2)




dev.off()
summary(lm(log(gamma_ratio)~log(Rec_Ad_ratio), data = ad_rec_gammas))
elastmod<-lm(log(gamma_ratio)~log(Rec_Ad_ratio), data = ad_rec_gammas)
confint(elastmod)


##Inner outer rather than islands####

inout_p<-fledge_p %>% 
  left_join(inout) %>% 
  group_by(in_out) %>% 
  mutate(sys_recap = mu_p + epsilonFY_p_fledge ) %>% 
  median_qi(sys_recap = invlogit(sys_recap), .width = 0.9) %>% 
  ungroup()

inout_p$sys_recap[inout_p$in_out==1]

outer_fledge<-stage_inout_summary %>% 
  filter(stage == "fledge" &in_out == 0)

inner_fledge<-stage_inout_summary %>% 
  filter(stage == "fledge" &in_out == 1)


corrected_fledge_outer<-outer_fledge$total/inout_p$sys_recap[inout_p$in_out==0]
corrected_fledge_outer
nest_surv_outer<-CMR_juv %>% 
  spread_draws(mu_phi) %>% 
  median_qi(nest_to_fledge_out =invlogit(mu_phi))

real_nestlings_outer<-corrected_fledge_outer/nest_surv_outer$nest_to_fledge_out



inner_fledge<-stage_inout_summary %>% 
  filter(stage == "fledge" &in_out == 1)

CMR_juv %>% 
  get_variables() %>% 
  head(n=20)
corrected_fledge_inner<-inner_fledge$total/inout_p$sys_recap[inout_p$in_out==1]
corrected_fledge_inner

nest_surv_inner<-CMR_juv %>% 
  spread_draws(mu_phi,phi_fledge_interac_inout) %>% 
  median_qi(nest_to_fledge_in =invlogit(mu_phi + phi_fledge_interac_inout))

real_nestlings_inner<-corrected_fledge_inner/nest_surv_inner$nest_to_fledge_in

real_nestlings_inner
real_nestlings_outer

adults_inout<-growth_factors %>% 
  group_by(in_out) %>% 
  summarise(tot_adults = sum(N))


adults_inout$real_nestlings <- c(real_nestlings_outer, real_nestlings_inner)
#Since we use N (yearly density) rather than unique adults to calculate birth rates
#we get the birth rate per year. This is because many individual are counted twice (they contributed to N in multiple years)
adults_inout<-adults_inout %>% 
  mutate(Birth_rate = real_nestlings/tot_adults) %>% 
  as.data.frame()
Birth_rate_outer<-adults_inout[1,4]
Fledge_rate_outer<-nest_surv_outer$nest_to_fledge_out

fledge_surv_outer<-CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,) %>% 
  median_qi(fledge_to_rec_out =invlogit(mu_phi +mu_rec_phi))

Rec_rate_outer<-fledge_surv_outer$fledge_to_rec_out

mean_ages<-cmr_raw %>% 
  group_by(in_out) %>% 
  mutate(Least_age = Least_age-1) %>% 
  mutate(Least_age_sq = Least_age^2) %>% 
  summarise(age = mean(Least_age)) %>% 
  ungroup()

inner_age<-mean_ages[2,2]$age
outer_age<-mean_ages[1,2]$age

ad_surv_outer<-CMR_adult_survival %>% 
  spread_draws(mu_phi, B_age, B_sex, B_age_sq,interac_age_sex) %>% 
  mutate(phi_outer = mu_phi + B_age*outer_age + B_age_sq*outer_age^2 +B_sex*0.5 +interac_age_sex*outer_age*0.5) %>% 
  median_qi(outer_ad_surv=invlogit(phi_outer))

Ad_rate_outer<-ad_surv_outer$outer_ad_surv

#Calculate elasticities
outer <- matrix(c(0, 0, Birth_rate_outer, Fledge_rate_outer, 0, 0, 0, Rec_rate_outer, Ad_rate_outer), nrow=3, byrow=TRUE,
                dimnames=list(c("nestling","fledgling","adult"),
                              c("nestling","fledgling","adult")))

outer
lambda(outer) #lambda is now 1 with the backcalculated nestlings! This better matches observed geometric lambda

#outer res from below. Calculate geometric mean.
prod(outer_res$lambda_obs, na.rm = T)^(1/(nrow(outer_res)-1)) #-1 because last row is na.

w_o <- eigen(outer)$vectors
v_o <- Conj(solve(w_o))
sensitivity(outer, zero = T)
(emat_out<-elasticity(outer))

#inner
Birth_rate_inner<-adults_inout[2,4]
Fledge_rate_inner<-nest_surv_inner$nest_to_fledge_in

fledge_surv_inner<-CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi, phi_rec_interac_inout) %>% 
  median_qi(fledge_to_rec_in =invlogit(mu_phi +mu_rec_phi + phi_rec_interac_inout))

Rec_rate_inner<-fledge_surv_inner$fledge_to_rec_in

ad_surv_inner<-CMR_adult_survival %>% 
  spread_draws(mu_phi, B_in_out, B_age, B_sex, B_age_sq,interac_age_sex) %>% 
  mutate(phi_inner = mu_phi + B_in_out + B_age*inner_age + B_age_sq*inner_age^2 +B_sex*0.5 +interac_age_sex*inner_age*0.5) %>% 
  median_qi(inner_ad_surv=invlogit(phi_inner))

Ad_rate_inner<-ad_surv_inner$inner_ad_surv

inner <- matrix(c(0, 0, Birth_rate_inner, Fledge_rate_inner, 0, 0, 0, Rec_rate_inner, Ad_rate_inner), nrow=3, byrow=TRUE,
                dimnames=list(c("nestling","fledgling","adult"),
                              c("nestling","fledgling","adult")))

inner
lambda(inner) #lambda is now 1 with the backcalculated nestlings! This matches the geometric mean of observed lambda

#inner res from below. Geometric mean.
prod(inner_res$lambda_obs, na.rm = T)^(1/(nrow(inner_res)-1)) #-1 because last row is na.

w_i <- eigen(inner)$vectors
v_i <- Conj(solve(w_i))


Re(eigen(inner)$values) #The built-in eigen function returns eigenvalues in descreasing order of magnitude or modulus.
sensitivity(inner, zero = T)
(emat_in<-elasticity(inner))

eig <- eigen(inner)
lambda <- Re(eig$values[1])  # Dominant eigenvalue

##Elasticity table####
in_out_list <- list()
in_out_list[[1]]<-emat_out
in_out_list[[2]]<-emat_in
elasticity_list_all<-c(in_out_list,elasticity_list)
names(elasticity_list_all)[[1]] <- "outer"
names(elasticity_list_all)[[2]] <- "inner"

simple_elast<-function(x){
  y<-data.frame(Vital_rate = c("$F$", "$\\phi_{ad}$"),
                Elasticity = c(x[1,3],x[3,3]))
  y
}
sublist<-lapply(elasticity_list_all, simple_elast)

elast_frame<-as.data.frame(do.call(rbind, sublist))

elast_frame$System_Island <- as.vector(rbind(names(elasticity_list_all), NA))

rownames(elast_frame) <- NULL

elast_frame<-elast_frame %>% 
  relocate(System_Island, .before = Vital_rate)
colnames(elast_frame) <- c("System/Island", "Vital rate", "Elasticity")
elast_frame

library(xtable)
output_filename <- "elasticity_table.tex"

xtable_result <- xtable(elast_frame, 
                        escape = FALSE, 
                        booktabs = TRUE)

print(
  xtable_result, 
  file = paste0("Tables empirical data/Tables output/", output_filename), 
  include.rownames = FALSE, 
  sanitize.text.function = identity, 
  sanitize.colnames.function = function(x) { paste0("\\multicolumn{1}{l}{", x, "}") },
  floating = FALSE
)
