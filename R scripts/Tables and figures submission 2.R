setwd("C:/Users/torhf/OneDrive - NTNU/Density regulation/R/Density regulation")

library(Matrix)
library(tidyverse)
library(tables)
library(knitr)
library(kableExtra)
library(readxl)
library(tidybayes)
library(arm)
library(xtable)

source("Functions/Functions to generate tables and figures from stan and brms outputs.R")
source("Functions/brm2table_short.r")
source("Functions/CMR2table_adults.r")
source("Functions/CMR2table_juveniles.R")
data_observed<-read.table("Cleaned data for density dependence analyses/cleaned_density_dependence_data.txt", sep=";", header=TRUE)



##Load latex variable names ####

paramlatex<-read_xlsx(path="Tables empirical data/Parameter to latex.xlsx")

##FIGURES####

###### Island slopes mean centred ####
#zmc_N (zNR)

density_variable = "zmc_N"
fitness = "r2s"
divide_intercepts = "yes"
average_metapars ="yes"
mean_centred = "yes"
r2s_zNR_new<-readRDS(file("Workspace backup/r2s_zNR_new.rds"))

summary(r2s_zNR_new)
brm2table_short(r2s_zNR_new, type="median", width = 0.9)

post_nb_full_Nm<-brm_blups_islands(r2s_zNR_new, density_variable, fitness, divide_intercepts,average_metapars)
post_nb_full_Nm


source("Functions/brm2table_short.r")

tab<-brm2table_short(r2s_zNR_new, type = "median", width = 0.9)

round(tab[tab$Parameter == "Sigma2_Isl_r0",c("Estimate",".lower", ".upper")],2)
var_slopes2<-round(tab[tab$Parameter == "Sigma2_Isl_gamma",c("Estimate",".lower", ".upper")],3)

var_slopes<-paste(var_slopes2[1], " (",var_slopes2[2],",", var_slopes2[3], ")", sep ="")

post_frame= post_nb_full_Nm

popsize_frame <- data_observed %>% 
  dplyr::select(Location, zNR_new) %>% 
  rename(popsize = zNR_new) %>% 
  rename(Island = Location)



distribution = "nb"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5
adjust_var = 1


source("Functions/Functions to generate tables and figures from stan and brms outputs.R")
island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

#Save v on n and w on N island plots to jpeg
jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100) 
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5)
# 2. Create a plot

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()

##### log N ####
density_variable = "n"
fitness = "r2s"
divide_intercepts= "yes"
average_metapars ="yes"
source("Functions/Functions to generate tables and figures from stan and brms outputs.R")
r2s_n_new<-readRDS(file("Workspace backup/r2s_n_new.rds"))

tab_n<-brm2table_short(r2s_n_new, type ="median" , width =0.9)

post_nb_full<-brm_blups_islands(r2s_n_new, density_variable, fitness, divide_intercepts,average_metapars)

min(exp(post_nb_full$r0))
max(exp(post_nb_full$r0))


var_slopes2<-round(tab_n[tab_n$Parameter == "Sigma2_Isl_gamma",c("Estimate",".lower", ".upper")],3)

var_slopes<-paste(var_slopes2[1], " (",var_slopes2[2],",", var_slopes2[3], ")", sep ="")
post_frame= post_nb_full

data_observed$n_new <- log(data_observed$N_corr)

popsize_frame <- data_observed %>% 
  dplyr::select(Location, n_new) %>% 
  rename(popsize = n_new) %>% 
  rename(Island = Location)
mean_centred = "no"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5
adjust_var = 1

source("Functions/Functions to generate tables and figures from stan and brms outputs.R")
island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)


distribution = "nb"
jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100)
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5)
# 2. Create a plot

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()


## Islands modelled as separate populations ####

#r2s zNR

r2s_zNR_single_island_mods_NEW<-readRDS(file("Workspace backup/r2s_zNR_single_island_mods_NEW.rds"))
r2s_zNR_single_island_mods_NEW[[1]]
divide_intercepts = "yes"

post_frame_single= get_single_island_estimates(r2s_zNR_single_island_mods_NEW,divide_intercepts, summary = TRUE)


var_slopes <- round(var(post_frame_single$g),3)

popsize_frame <- data_observed %>% 
  dplyr::select(Location, zNR_new) %>% 
  rename(popsize = zNR_new) %>% 
  rename(Island = Location)


mean_centred ="yes"

density_variable ="zmc_N"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5
fitness = "r2s"
adjust_var = 1
island_slopes(post_frame_single,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

distribution = "nb"

jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_single_islands_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100) 
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_single_islands_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5) 
# 2. Create a plot

island_slopes(post_frame_single,popsize_frame,data_observed,mean_centred, density_variable,var_slopes,adjust_var, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()


##TABLES####
######r2S, recruits and adult survival CMR N_m table####
r2s_zNR_new<-readRDS(file("Workspace backup/r2s_zNR_new.rds"))
N_recruits_zNR_new<-readRDS(file("Workspace backup/N_recruits_zNR_new.rds"))
CMR_nest_adult_survival<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")
round(summary(CMR_nest_adult_survival)$summary[c(1:21),c(6, 4, 8)], 3)

CMR_nest_adult_survival %>% 
  spread_draws(mu_p) %>% 
  median_qi(mean_p = invlogit(mu_p), .width =0.9)

width =0.9
r2s_zNR_new_table <-brm2table_short(r2s_zNR_new, type ="median", width)
N_recruits_zNR_new_table<-brm2table_short(N_recruits_zNR_new, type = "median", width)


CMR_nest_adult_posterior<-as.data.frame(CMR_nest_adult_survival)

CMR_nest_adult_posterior<-CMR_nest_adult_posterior %>% 
  dplyr::select(any_of(paramlatex$Parameter))

CMR_nest_adult_posterior<-CMR_nest_adult_posterior %>% 
  mutate(mu_phi_rec_pred = mu_phi + mu_phi_rec) %>% 
  mutate(B_in_out_rec_pred = mu_phi_rec_pred + B_in_out_rec) %>% 
  mutate(B_in_out_rec_diff = mu_phi_rec_pred - B_in_out_rec_pred) %>% 
  dplyr::select(-mu_phi_rec, -B_in_out_rec, -B_in_out_rec_pred) %>% 
  rename(mu_phi_rec = mu_phi_rec_pred) %>% 
  rename(B_in_out_rec = B_in_out_rec_diff)

CMR_nest_adult_posterior_long<-CMR_nest_adult_posterior %>% 
  pivot_longer(1:ncol(CMR_nest_adult_posterior), names_to = "Parameter", values_to = "Estimate")



CMR_nest_adult_table<-CMR_nest_adult_posterior_long %>% 
  group_by(Parameter) %>% 
  median_qi(Estimate, .width = width)

CMR_nest_adult_table

CMR_nest_table<-CMR_nest_adult_table %>% 
  filter(str_detect(Parameter , "rec"))

CMR_adult_table<-CMR_nest_adult_table %>% 
  filter(str_detect(Parameter , "rec", negate = T))

roundinvlogit<-function(x){round(invlogit(x), 2)
}
mu_p<-CMR_nest_adult_table %>% 
  filter(Parameter == "mu_p") %>% 
  mutate(across(where(is.numeric), roundinvlogit))

Sigma2_YI_p<-CMR_nest_adult_table %>% 
  filter(Parameter == "Sigma2_YI_p") %>% 
  mutate(across(where(is.numeric), ~round(.x,2)))


table_capt<-"The effects of density (mean-centred and scaled adult density) and 
other sources of variation on house sparrow individual demographic contributions 
and its components; the number of recruits produced, nestling to adult survival and adult survival, with
90\\% credible intervals in parentheses."

tbl_label <- "Nm_nb_tbl"
table_name <- "results_r2s_rec_nest_ad_adult"
table_footnote <- paste0("The average recapture probability was p = ",mu_p$Estimate, " (90\\% CI = ", mu_p$.lower, ", ", mu_p$.upper, ") with a variance \\\ among years within populations on the latent scale of ", Sigma2_YI_p$Estimate, " (90\\% CI = ", Sigma2_YI_p$.lower, ", ", Sigma2_YI_p$.upper, "). Island, Island:density 
and Island cov represent the (co)variances of island-specific intercepts and slopes 
for the effect of density. Individual and Island-year represent variances in individual and Island-year 
specific intercepts. The unexplained variance in adult survival is not presented since 
it is captured by the Bernoulli process.")


tables <- list()

tables[[1]] <- r2s_zNR_new_table
tables[[2]] <- N_recruits_zNR_new_table
tables[[3]] <- CMR_nest_table
tables[[4]] <- CMR_adult_table

col_headers <- c(
  "\\textbf{Demo. contribution}",
  "\\textbf{Recruits}",
  "\\textbf{Nest.to ad. surv.}",
  "\\textbf{Ad. surv.}"
)
sigfigs = 2
n_fixef = 8
fit_table<-table_wrangler_smart(tables,col_headers,sigfigs,n_fixef, paramlatex)
fit_table<-fit_table %>% 
  filter(str_detect( `\\textbf{Parameter}`,"Mean recapture probability|Sigma2_YI_p", negate = T))
fit_table

#results_frame = fit_table
footn <- list(pos=list(0), command= NULL)
footn$pos[[1]] <- c(nrow(fit_table))
footn$command <- c(paste0("\\hline\\multicolumn{4}{l}{Footnote: See Appendix S5 for recapture probabilities. Habitat and Density:habitat represent}\\\\
   \\multicolumn{4}{l}{the effect of habitat system and its interaction with density. Age:Sex represents the interaction}\\\\
   \\multicolumn{4}{l}{of age and sex. Island, Island:density and Island cov represent the (co)variances of island-}\\\\
                          \\multicolumn{4}{l}{specific intercepts and slopes for the effect of density. Individual and Island-year represent}\\\\
                          \\multicolumn{4}{l}{variances in individual and Island-year specific intercepts. The unexplained variance in adult}\\\\
                          \\multicolumn{4}{l}{survival is not presented since it is captured by the Bernoulli process.}", "\n", sep = ""))


#"\\hline","\\multicolumn{4}{l}{Footnote: The average recapture probability was p = ", mu_p$Estimate, "(90\\% CI = ", mu_p$.lower, ", ", mu_p$.upper, " with a variance}", "\\\\",
#"\\multicolumn{4}{l}{among years within populations on the latent scale of ", Sigma2_YI_p$Estimate,  "(90\\% CI = ", Sigma2_YI_p$.lower, "," , Sigma2_YI_p$.upper,"). Island,}", "\\\\",
#"\\multicolumn{4}{l}{Island:density and Island cov represent the (co)variances of island-specific intercepts and}", "\\\\",
#"\\multicolumn{4}{l}{slopes for the effect of density. Individual and Island-year represent variances in individual}", "\\\\",
#"\\multicolumn{4}{l}{and Island-year specific intercepts. The unexplained variance in adult survival is not}", "\\\\",
#"\\multicolumn{4}{l}{presented since it is captured by the Bernoulli process.}"



#add.to.row = footn, hline.after = c(-1,0),
print(xtable(fit_table,type = "latex", caption = table_capt, label = tbl_label,floating=FALSE, escape = F,booktabs = T), file = paste0("Tables empirical data/Tables output/",table_name,".tex"),caption.placement = "top", include.rownames=FALSE,
      add.to.row = footn, hline.after = c(-1,0),
      sanitize.text.function = identity, sanitize.colnames.function=function(x){paste0("\\multicolumn{1}{l}{",x,"}")},size="\\fontsize{12pt}{20pt}\\selectfont") #multicolumn{1}{l} the l left-aligns the col headers



# To make it simpler we create captions in footnotes in Overleaf

xtable_result <- xtable(fit_table, 
                        floating = FALSE, 
                        escape = FALSE, 
                        booktabs = TRUE)

print(
  xtable_result,
  file = paste0("Tables empirical data/Tables output/", table_name,".tex"),
  include.rownames = FALSE,
  floating = FALSE,
  sanitize.text.function = identity,
  sanitize.colnames.function = function(x) {
    paste0("\\multicolumn{1}{l}{", x, "}")
  }
)


##Descriptive statistics in chronological order ####

#####Mean popsizes
demo_N_corr <- data_observed %>% 
  dplyr::select(Location, Year, N_corr) %>% 
  unique()

mean(demo_N_corr$N_corr)

demo_N_corr %>% 
  group_by(Location) %>% 
  summarise(mean(N_corr))

#####Change in population size from year to year, proportions and individuals ####

popsizes<-read.table("Cleaned data for density dependence analyses/estimated_pop_sizes_nestling_to_adult_model_with_lurøy_onøy.txt",sep =";", encoding = "UTF-8", header = T)

popsizes %>% 
  group_by(Location) %>% 
  summarise(min(Year))

popsizes %>% 
  filter(Location == "aldra")

head(popsizes)
isl_pop<-popsizes %>% 
  dplyr::select(Location, Year, N_corr) %>% 
  unique() %>% 
  arrange(Location, Year)


isl_pop %>% 
  mean_qi(grandmean_popsize = N_corr, .width = 0.9)
isl_pop %>% 
  group_by(Location) %>% 
  mean_qi(mean_popsize = N_corr, .width = 0.9)

min(isl_pop$N_corr)
max(isl_pop$N_corr)
isl_pop_prop<-isl_pop %>% 
  group_by(Location) %>% 
  mutate(N_plusone = lead(N_corr)) %>% 
  filter(Year<2021) %>% 
  filter(N_corr != 0) %>% 
  mutate(deltaN = N_plusone/N_corr) %>% 
  mutate(N_change = abs(1- deltaN)) %>%  #abs returns absolute values
  mutate(Mean_isl_change = mean(N_change)) %>% 
  dplyr::select(Location, Mean_isl_change) %>% 
  unique()
min(isl_pop_prop$Mean_isl_change)
max(isl_pop_prop$Mean_isl_change)

isl_pop_inds<-isl_pop %>% 
  group_by(Location) %>% 
  mutate(N_plusone = lead(N_corr)) %>% 
  filter(Year<2021) %>% 
  mutate(N_change_inds = N_plusone-N_corr) %>% 
  mutate(N_change_inds = abs(N_change_inds)) %>%  #abs returns absolute values
  mean_qi(Mean_N_change_inds = N_change_inds, .width = 0.9)

isl_pop_inds
min(isl_pop_inds$Mean_N_change_inds)
max(isl_pop_inds$Mean_N_change_inds)

isl_pop %>% 
  group_by(Location) %>% 
  mutate(N_plusone = lead(N_corr)) %>% 
  filter(Year<2021) %>% 
  mutate(N_change_inds = N_plusone-N_corr) %>% 
  mutate(N_change_inds = abs(N_change_inds)) %>%  #abs returns absolute values
  ungroup() %>% 
  mean_qi(Mean_N_change_inds = N_change_inds, .width = 0.9)

#####Results r2s model ####
get_variables(r2s_zNR_new)[1:12]

#mean strength of density dependence metapop
r2s_zNR_new %>% 
  spread_draws(b_zNR_new, `b_zNR_new:in_out`) %>% 
  mutate(zNR_metapop = b_zNR_new + `b_zNR_new:in_out`*0.5) %>% 
  median_qi(zNR_metapop, .width=0.9)

#reduction in fitness with increase of 1 zNR (40 inds) metapop

r2s_zNR_new  %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR_new, `b_zNR_new:in_out`) %>%
  mutate(meta_g = b_zNR_new + `b_zNR_new:in_out`*0.5) %>% 
  mutate(fitness_equil = b_Intercept + b_in_out*0.5) %>% 
  mutate(fitness_1zNR = fitness_equil + meta_g*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(1-prop, .width = 0.9)

#reduction in fitness outer system increase of 1 zNR

r2s_zNR_new %>% 
  spread_draws(b_Intercept, b_zNR_new) %>%
  mutate(fitness_equil = b_Intercept) %>% 
  mutate(fitness_1zNR = fitness_equil + b_zNR_new*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(1-prop, .width = 0.9)

#reduction in fitness inner system increase of 1 zNR
r2s_zNR_new %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR_new, `b_zNR_new:in_out`) %>%
  mutate(inner_g = b_zNR_new + `b_zNR_new:in_out`) %>% 
  mutate(fitness_equil = b_Intercept + b_in_out) %>% 
  mutate(fitness_1zNR = fitness_equil + inner_g*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(1-prop, .width = 0.9)
#variance in slopes
r2s_zNR_new %>% 
  spread_draws(sd_Location__zNR_new) %>% 
  median_qi(sd_Location__zNR_new^2, .width = 0.9)

#mean r2s outer

r2s_zNR_new %>% 
  spread_draws(b_Intercept) %>%
  mutate(fitness_equil = exp(b_Intercept-log(2))) %>% 
  median_qi(fitness_equil,.width = 0.9)

#mean r2s inner
r2s_zNR_new %>% 
  spread_draws(b_Intercept,b_in_out) %>%
  mutate(fitness_equil = exp(b_Intercept + b_in_out-log(2))) %>% 
  median_qi(fitness_equil, .width = 0.9)

#####adult survival ####
CMR_nest_adult_survival<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")
get_variables(CMR_nest_adult_survival)[1:20]

#survival
CMR_nest_adult_survival %>% 
  spread_draws(mu_phi, B_in_out, B_sex)%>%
  mutate(phi_inner = invlogit(mu_phi + B_in_out + B_sex*0.5)) %>% 
  mutate(phi_outer = invlogit(mu_phi + B_sex*0.5)) %>% 
  median_qi(phi_outer,phi_inner, .width =0.9)

#inner gamma
CMR_nest_adult_survival %>% 
  spread_draws(gamma, gamma_in_out) %>%
  mutate(gamma_inner = gamma + gamma_in_out) %>% 
  median_qi(gamma,gamma_inner, .width = 0.9)

##Difference survival mean aged adult inner and outer
CMR_nest_adult_survival %>% 
  spread_draws(B_in_out) %>% 
  median_qi(B_in_out, .width =0.9)

CMR_nest_adult_survival %>% 
  spread_draws(mu_phi,B_sex, B_in_out) %>% 
  mutate(inner_average_phi = 
           invlogit(mu_phi + B_sex*0.5  + B_in_out)) %>% 
  mutate(outer_average_phi = 
           invlogit(mu_phi + B_sex*0.5)) %>% 
  median_qi(diff_surv  =inner_average_phi-outer_average_phi, .width = 0.9) %>% 
  mutate(diff_surv = round(diff_surv,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

#####nestling to adult survival####

#inner gamma
CMR_nest_adult_survival %>% 
  spread_draws(gamma_rec, gamma_in_out_rec) %>% 
  mutate(inner_gamma_rec = gamma_rec +gamma_in_out_rec) %>% 
  median_qi(inner_gamma_rec, .width =0.9)

#nestling to adult survival mean metapop
CMR_nest_adult_survival %>% 
  spread_draws(mu_phi,mu_phi_rec, B_in_out_rec) %>% 
  mutate(mu_phi_rec_tot = mu_phi +mu_phi_rec + B_in_out_rec*0.5) %>% 
  median_qi(invlogit(mu_phi_rec_tot), .width =0.9)


#nestling to adult survival mean metapop
CMR_nest_adult_survival %>% 
  spread_draws(mu_phi,mu_phi_rec, B_in_out_rec) %>% 
  mutate(mu_phi_rec_outer = invlogit(mu_phi +mu_phi_rec)) %>% 
  mutate(mu_phi_rec_inner = invlogit(mu_phi +mu_phi_rec + B_in_out_rec)) %>% 
  median_qi(mu_phi_rec_outer, mu_phi_rec_inner, .width =0.9)

######nestling to fledge, fledge to rec ####

CMR_juveniles<-readRDS("Workspace backup/CMR_juvenile_survival_zNR_newpop.rds")


##The file is very big. Extract the parameter of intereset

CMR_juv_sub<-CMR_juveniles%>% 
  spread_draws(mu_phi, phi_fledge_interac_inout, mu_rec_phi, phi_rec_interac_inout,gamma,gamma_fledge_interac_inout,  gamma_rec, gamma_rec_interac_inout)


##gamma diff
CMR_juv_sub %>% 
  spread_draws(gamma_rec_interac_inout) %>% 
  median_qi(gamma_rec_interac_inout, .width = 0.9)
CMR_juv_sub %>% 
  spread_draws(gamma_fledge_interac_inout) %>% 
  median_qi(gamma_fledge_interac_inout, .width = 0.9)

##differences inner outer phi

CMR_juv_sub %>% 
  median_qi(phi_fledge_interac_inout, phi_rec_interac_inout, .width =0.9)

##phi fledge survival prop diff
CMR_juv_sub %>% 
  spread_draws(mu_phi, phi_fledge_interac_inout) %>% 
  mutate(inner_fledge_intercept = mu_phi + phi_fledge_interac_inout) %>% 
  mutate(prop_fledge_phi = invlogit(mu_phi)/invlogit(inner_fledge_intercept)) %>% 
  median_qi(prop_fledge_phi, .width =0.9)

##phi fledge percentate point survival diff
CMR_juv_sub %>% 
  spread_draws(mu_phi, phi_fledge_interac_inout) %>% 
  mutate(inner_fledge_intercept = mu_phi + phi_fledge_interac_inout) %>% 
  mutate(diff_fledge_phi_percentpoints = invlogit(mu_phi) - invlogit(inner_fledge_intercept)) %>% 
  median_qi(diff_fledge_phi_percentpoints, .width =0.9)


##predictions inner and outer phi

CMR_juv_sub %>%
  mutate(inner_fledge_intercept = invlogit(mu_phi + phi_fledge_interac_inout)) %>% 
  mutate(outer_fledge_intercept = invlogit(mu_phi)) %>%
  median_qi(inner_fledge_intercept, outer_fledge_intercept, .width =0.9)

CMR_juv_sub %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi)) %>% 
  median_qi(inner_rec_phi, outer_rec_phi, .width =0.9)

##phi rec survival prop diff
CMR_juv_sub %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi)) %>%
  mutate(prop_rec_phi = inner_rec_phi/outer_rec_phi) %>% 
  median_qi(prop_rec_phi, .width = 0.9)

CMR_juv_sub %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi)) %>%
  mutate(diff_rec_phi_percentpoints = inner_rec_phi - outer_rec_phi) %>% 
  median_qi(diff_rec_phi_percentpoints, .width = 0.9)


#survival to adulthood, contrast to nestling to adult model
CMR_juv_sub %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi)) %>% 
  mutate(outer_fledge_phi =  invlogit(mu_phi)) %>% 
  mutate(inner_fledge_phi =  invlogit(mu_phi + phi_fledge_interac_inout)) %>% 
  mutate(outer_nest_rec_phi = outer_fledge_phi * outer_rec_phi) %>% 
  mutate(inner_nest_rec_phi = inner_fledge_phi * inner_rec_phi) %>% 
  median_qi(outer_nest_rec_phi, inner_nest_rec_phi, .width =0.9)

#gamma fledge
CMR_juv_sub %>% 
  mutate(gamma_fledge_inner = gamma + gamma_fledge_interac_inout) %>% 
  median_qi(gamma,gamma_fledge_inner, .width = 0.9)


##gamma rec
head(CMR_juv_sub)
CMR_juv_sub %>% 
  median_qi(gamma_rec)

CMR_juv_sub %>% 
  mutate(gamma_rec_inner = gamma_rec + gamma_rec_interac_inout) %>% 
  median_qi(gamma_rec,gamma_rec_inner, .width = 0.9)

######Nestling production ####
nestlings_zNR_nb_brm_noindslope_new<-readRDS(file("Workspace backup/nestlings_zNR_nb_brm_noindslope_new.rds")) 
nestlings_zNR_nb_brm_noindslope_new %>% 
  spread_draws(b_Intercept) %>% 
  mutate(outer_nestprod = exp(b_Intercept -log(2))) %>% 
  median_qi(outer_nestprod, .width =0.9)

nestlings_zNR_nb_brm_noindslope_new %>% 
  spread_draws(b_Intercept, b_in_out) %>% 
  mutate(inner_nestprod = exp(b_Intercept+ b_in_out -log(2))) %>% 
  median_qi(inner_nestprod, .width =0.9)

nestlings_zNR_nb_brm_noindslope_new %>% 
  spread_draws(b_Intercept, b_in_out) %>% 
  mutate(meta_nestprod = exp(b_Intercept+ b_in_out*0.5 -log(2))) %>% 
  median_qi(meta_nestprod, .width =0.9)

nestlings_zNR_nb_brm_noindslope_new %>% 
  spread_draws(b_zNR_new) %>% 
  mutate(outer_gamma = b_zNR_new) %>% 
  median_qi(outer_gamma, .width =0.9)


nestlings_zNR_nb_brm_noindslope_new %>% 
  spread_draws(b_zNR_new, `b_zNR_new:in_out`) %>% 
  mutate(inner_gamma = b_zNR_new + `b_zNR_new:in_out`) %>% 
  median_qi(inner_gamma, .width =0.9)





##FIGURES SUPPORTING MATERIAL ####

####### Correlation between estimated pop size (pre-breeding census) and our fitness estimates (r2s/2) ####
fits<-data_observed %>% 
  group_by(Location, Year) %>% 
  mutate(tot_fitness = sum(r2s)/2) %>% 
  ungroup %>% 
  dplyr::select(Location, Year,tot_fitness) %>% 
  mutate(Year = Year +1 ) %>% 
  filter(Year != 2020) %>% 
  unique()

pops <- data_observed %>% 
  dplyr::select(Location, Year, N_corr)

pops <- data_observed %>% 
  dplyr::select(Location, Year, N_corr) %>% 
  #filter(Year != 1998) %>% 
  unique()
fits<-fits %>% 
  left_join(pops, by = c("Location", "Year")) %>% 
  arrange(Location, Year)

yrs<-seq(1999, 2019, by = 4)

pdf("Figures empirical data/estimated pop size vs RS.pdf",
width = 6.7,
height = 7.1)

par(mfrow=c(6,2))
par(mar = c(1,2,3,1))
par(oma = c(4,4,1,1))
for(i in 1:11){
  single_island <- subset(fits, Location == unique(fits$Location)[i])
  
  plot(N_corr~Year, data = single_island, main = unique(single_island$Location), xlab = "", ylab ="", type ="l", xaxt = "n", ylim = c(0, max(N_corr)))
  axis(1, at=yrs)
  points(tot_fitness~Year, data = single_island, col = "red", type = "l")
} 

plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab="", ylab="")
legend('center',
       legend = c("Estimated pop size", "Recruitment + survival"), 
       col = c("black","red"), 
       lwd=2, xpd = TRUE, cex = 2, seg.len=1, bty = 'n')
mtext(text="Year",side=1,line=2,outer=TRUE)
mtext(text="Number of adults",side=2,line=1,outer=TRUE)

dev.off()

##  Correlation between pop size and gamma and r0 ####
r2s_zNR_new<-readRDS(file("Workspace backup/r2s_zNR_new.rds"))
r2s_zNR_single_island_mods_NEW<-readRDS(file("Workspace backup/r2s_zNR_single_island_mods_NEW.rds"))
r2s_n_new<-readRDS(file("Workspace backup/r2s_n_new.rds"))

density_variable = "n"
fitness = "r2s"
divide_intercepts= "yes"
average_metapars ="yes"
source("Functions/Functions to generate tables and figures from stan and brms outputs.R")

post_nb_full<-brm_blups_islands(r2s_n_new, density_variable, fitness, divide_intercepts,average_metapars)

post_nometa_n<-post_nb_full %>% 
  filter(Island != "metapopulation") %>% 
  rename(Location = Island)

island_N_means<-data_observed %>% 
  dplyr::select(Location,in_out,Mean_N_corr_isl) %>% 
  unique() 
get_variables(r2s_zNR_new)[1:10]

int_slope<-r2s_zNR_new %>% 
  spread_draws(b_Intercept, b_in_out, b_zNR_new,`b_zNR_new:in_out`,
               r_Location[Location, variable])

int_slope$Location <- gsub("indre.kvarøy", "indre kvarøy", int_slope$Location)
unique(int_slope$Location)




int_slope<-int_slope %>% 
  left_join(island_N_means)

int_slope<-int_slope %>% 
  pivot_wider(names_from = "variable", values_from = r_Location)

int_slope<-int_slope %>% 
  mutate(r_mu = ifelse(in_out == 0, b_Intercept + Intercept -log(2), b_Intercept +b_in_out + Intercept -log(2))) %>%
  mutate(gamma = ifelse(in_out == 0, b_zNR_new + zNR_new, b_zNR_new +`b_zNR_new:in_out` + zNR_new))

r0_backtransformed_zNR<-int_slope %>% 
  mutate(r0 = r_mu - (Mean_N_corr_isl/40)*gamma) %>% 
  group_by(Location) %>% 
  median_qi(r0,gamma, .width = 0.9)

max(exp(r0_backtransformed_zNR$r0))

r0_backtransformed_zNR %>% 
  arrange(gamma)

stopifnot(
identical(island_N_means$Location, names(r2s_zNR_single_island_mods_NEW))
)
res_list<-list()
for(i in 1:11){
single_mod<-r2s_zNR_single_island_mods_NEW[[i]]

single_pop_N<-island_N_means[i,]

single_mod<-single_mod %>% 
  spread_draws(b_Intercept, b_gamma) %>% 
  mutate(r0 = b_Intercept- b_gamma*single_pop_N$Mean_N_corr_isl/40 -log(2)) %>% 
  median_qi(r0, b_gamma)
single_mod$Location <-single_pop_N$Location
res_list[[i]]<-single_mod
res_list
}
zNR_single_backtransformed<-do.call(rbind.data.frame, res_list)

r0_backtransformed_zNR
zNR_single_backtransformed
post_nometa_n

gr0_zNr_meta<-r0_backtransformed_zNR %>% 
  left_join(island_N_means)


gr0_zNr_single<-zNR_single_backtransformed %>% 
  left_join(island_N_means)

gr0n_meta<-post_nometa_n %>% 
  left_join(island_N_means)



pdf("Figures empirical data/gamma r0 mean densities.pdf",
    width = 10,
    height = 6)
par(mfrow = c(3,3))
par(mar= c(2,2,2,2))
par(oma= c(6,1,4,2))

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression('Logistic, meta-population model (N'[m]*')'), cex = 1.4)

plot(gamma ~Mean_N_corr_isl, data = gr0_zNr_meta, xlab="", ylab="",cex = 1.5,
     col=as.factor(gr0_zNr_meta$in_out),pch= as.numeric(as.factor(gr0_zNr_meta$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0_zNr_meta$in_out)), col = 1:2, cex = 0.1, pch=unique(as.numeric(as.factor(gr0_zNr_meta$in_out))))
plot(r0~Mean_N_corr_isl, data = gr0_zNr_meta, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_meta$in_out),pch= as.numeric(as.factor(gr0_zNr_meta$in_out)), cex.axis=1.3)

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression('Logistic, single population model (N'[m]*')'), cex = 1.4)

plot(b_gamma~Mean_N_corr_isl, data = gr0_zNr_single, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_single$in_out), pch= as.numeric(as.factor(gr0_zNr_single$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0_zNr_single$in_out)), col = 1:2, cex = 1, pch=unique(as.numeric(as.factor(gr0_zNr_single$in_out))))
plot(r0~Mean_N_corr_isl, data = gr0_zNr_single, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_single$in_out), pch= as.numeric(as.factor(gr0_zNr_single$in_out)), cex.axis=1.3)

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression("Gompertz meta-population model (n)"), cex = 1.4)


plot(g~Mean_N_corr_isl, data = gr0n_meta, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0n_meta$in_out),pch= as.numeric(as.factor(gr0n_meta$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0n_meta$in_out)),  col = 1:2, cex = 0.1, pch=unique(as.numeric(as.factor(gr0n_meta$in_out))))
plot(r0~Mean_N_corr_isl, data = gr0n_meta, xlab="", ylab="",cex = 1.5,
     col=as.factor(gr0n_meta$in_out),pch= as.numeric(as.factor(gr0n_meta$in_out)), cex.axis=1.3)
legend("bottomright",
       legend = c("Inner farm","Outer non-farm"), 
       col = c("red", "black"), 
       xpd = TRUE, cex = 1.2,pch=unique(as.numeric(as.factor(gr0n_meta$in_out))), seg.len=1)

mtext(expression(gamma), side = 3, outer = T, cex = 1.3)
mtext(expression("r"[0]), side = 3, outer = T, adj = 0.85, cex = 1.3)
mtext("Island-specific mean density", side = 1, outer = T, adj = 0.7, cex = 1.3,line = 2)

dev.off()


##TABLES SUPPORTING MATERIAL####

## Juvenile survival CMR, nestling production N_m ####

CMR_juveniles<-readRDS("Workspace backup/CMR_juvenile_survival_zNR_newpop.rds")
nestlings_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/nestlings_zNR_nb_brm_noindslope_new.rds")) 
summary(nestlings_zNR_nb_brm_noindslope)

width = 0.9
nestling_prod_res <- brm2table_short(nestlings_zNR_nb_brm_noindslope, type = "median",width)

reverse_inner_outer="no"
add_difference_to_intercept="yes"

juv_res<-CMR2table_juveniles(CMR_juveniles,paramlatex,reverse_inner_outer,add_difference_to_intercept, width)


(nest2fledge_res<-juv_res[[1]])
(fledge2rec_res<-juv_res[[2]])
nestling_prod_res


roundinvlogit<-function(x){round(invlogit(x), 2)
}
mu_p_fledge<-nest2fledge_res %>% 
  filter(Parameter == "mu_p_fledge") %>% 
  mutate(across(where(is.numeric), roundinvlogit))

Sigma2_YI_p_fledge<-nest2fledge_res %>% 
  filter(Parameter == "Sigma2_YI_p_fledge") %>% 
  mutate(across(where(is.numeric), ~round(.x,2)))

mu_rec_p<-fledge2rec_res %>% 
  filter(Parameter == "mu_rec_p") %>% 
  mutate(across(where(is.numeric), roundinvlogit))

Sigma2_YI_p_rec<-fledge2rec_res %>% 
  filter(Parameter == "Sigma2_YI_p_rec") %>% 
  mutate(across(where(is.numeric), ~round(.x,3)))


raw_tables<-list(nestling_prod_res,nest2fledge_res, fledge2rec_res)


col_headers<-c('\\textbf{No. of nestlings produced}', '\\textbf{Nestling survival}', '\\textbf{Fledgling survival}')

n_fixef = 8
sigfigs = 2
fit_table<-table_wrangler_smart(raw_tables,col_headers,sigfigs,n_fixef, paramlatex)
fit_table<-fit_table %>% 
  filter(str_detect(`\\textbf{Parameter}`,"Mean recapture probability|Sigma2_YI_p_fledge|Sigma2_YI_p_rec", negate = T))
fit_table
#results_frame$`Nestling to fledgling`[results_frame$`Nestling to fledgling` =="NA (NA, NA)"] <- NA
#results_frame$`Fledgling to recruit`[results_frame$`Fledgling to recruit` =="NA (NA, NA)"] <- NA  

table_capt <- "Medians for posterior distributions for the effects of density
(mean-centred and scaled adult density) and other sources of variation 
on number of nestlings produced by successfully reproducing parents (negative binomial), 
nestling to fledgling survival and fledgling to recruit survival (Bernoulli, 
capture-mark-recapture) with 90\\% credible intervals in parentheses."

tbl_label <- "surv_nest_fledge_rec_rep"
table_name <- "surv_nest_fledge_rec_rep"
table_footnote <- paste0("Nestling to fledgling and fledgling to recruit survival was modelled 
                         in a mark-recapture framework. The average fledgling recapture 
                         probability was p = ",mu_p_fledge$Estimate, " 
                         (90\\% CI = ", mu_p_fledge$.lower, ", ", mu_p_fledge$.upper, ") 
                         with a variance among years within populations on the latent scale of 
                         ", Sigma2_YI_p_fledge$Estimate, " 
                         (90\\% CI = ", Sigma2_YI_p_fledge$.lower, ", ", Sigma2_YI_p_fledge$.upper, ").",
                         "The average recruit recapture 
                         probability was p = ",mu_rec_p$Estimate, " 
                         (90\\% CI = ", mu_rec_p$.lower, ", ", mu_rec_p$.upper, ") 
                         with a variance among years within populations on the latent scale of 
                         ", Sigma2_YI_p_rec$Estimate, " 
                         (90\\% CI = ", Sigma2_YI_p_rec$.lower, ", ", Sigma2_YI_p_rec$.upper, ").
                         Juvenile survival was estimated in a single model but here full parameter 
estimates of fledgling survival are presented rather than the difference in survival
between nestlings and fledglings. Island, Island:density and Island cov represent
the (co)variances of island-specific intercepts and slopes for the effect of density. Individual and Island-year 
represent variances in individual and Island-year specific intercepts.
The unexplained variance in Bernoulli models is not presented since 
it is captured by the Bernoulli process.")



footn <- list(pos=list(0), command= NULL)
footn$pos[[1]] <- c(nrow(fit_table))
footn$command <- c(paste0("\\hline","\\multicolumn{4}{l}{Footnote: See Appendix S5 for recapture probabilities. Island, Island:density and Island cov}","\\\\",
                          "\\multicolumn{4}{l}{represent the (co)variances of island-specific intercepts and slopes for the effect of density. Individual}","\\\\",
                          "\\multicolumn{4}{l}{and Island-year represent variances in individual and Island-year specific intercepts. The unexplained }","\\\\",
                          "\\multicolumn{4}{l}{variance in Bernoulli models is not presented since it is captured by the Bernoulli process.}"))

#"\\hline","\\multicolumn{4}{l}{Footnote: Nestling and fledgling survival was modelled in a mark-recapture framework. The}", "\\\\", 
#"\\multicolumn{4}{l}{average fledgling recapture was p = ", mu_p_fledge$Estimate,  " (90\\% CI = ", mu_p_fledge$.lower, ", ", mu_p_fledge$.upper, ") with a variance among years}", "\\\\", 
#"\\multicolumn{4}{l}{within populations on the latent scale of ", Sigma2_YI_p_fledge$Estimate, " (90\\% CI = ", Sigma2_YI_p_fledge$.lower, ", ", Sigma2_YI_p_fledge$.upper, "). The average recruit recapture}","\\\\",
#"\\multicolumn{4}{l}{probability was p = ",mu_rec_p$Estimate, " (90\\% CI = ", mu_rec_p$.lower, ", ", mu_rec_p$.upper, ") with a variance among years within populations}", "\\\\",
#"\\multicolumn{4}{l}{on the latent scale of ", Sigma2_YI_p_rec$Estimate," (90\\% CI = ", Sigma2_YI_p_rec$.lower, ", ", Sigma2_YI_p_rec$.upper, "). Island, Island:density and Island cov represent}", "\\\\", 
#"\\multicolumn{4}{l}{the (co)variances of island-specific intercepts and slopes for the effect of density. Individual and}", "\\\\", 
#"\\multicolumn{4}{l}{Island-year represent variances in individual and Island-year specific intercepts.The unexplained}", "\\\\",
#"\\multicolumn{4}{l}{variance in Bernoulli models is not presented since it is captured by the Bernoulli process.}"))


#add.to.row = footn, hline.after = c(-1,0),
print(xtable(fit_table,type = "latex", caption = table_capt, label = tbl_label,floating=FALSE, escape = F,booktabs = T), file = paste0("Tables empirical data/Tables output/",table_name,".tex"),caption.placement = "top", include.rownames=FALSE,
      add.to.row = footn, hline.after = c(-1,0),
      sanitize.text.function = identity, sanitize.colnames.function=function(x){paste0("\\multicolumn{1}{l}{",x,"}")},size="\\fontsize{12pt}{20pt}\\selectfont") #multicolumn{1}{l} the l left-aligns the col headers


##r2s log N table####
r2s_n_new<-readRDS(file("Workspace backup/r2s_n_new.rds"))
n_table<-brm2table_short(r2s_n_new, type ="median" , width =0.9)

table_capt<-"Medians for posterior distributions for the effects of log density (log adult density) and other sources of
variation on individual demographic contribution (recruit production + adult survival, with 90% credible intervals
in parentheses. Habitat and Density:habitat represent the effect of habitat system and its interaction with density.
Age:Sex represents the interaction of age and sex. Island, Island:density and Island cov represent the (co)variances
of island- specific intercepts and slopes for the effect of density. Individual and Island-year represent variances in
individual and Island-year specific intercepts. The Shape parameter captures the overdispersion in negative binomial
models."

tables <- list()

tables[[1]] <- n_table

col_headers<-c("\\textbf{Demographic contribution}")
tbl_label <- "n_nb_tbl"
table_name <- "results_fixed_n_noindslope"
table_footnote <-NULL
sigfigs = 2
n_fixef = 8
fit_table<-table_wrangler_smart(tables,col_headers,sigfigs,n_fixef, paramlatex)

res_tables(fit_table, table_capt, tbl_label,table_footnote, table_name, n_fixef)

popsizes<-data_observed %>% 
  dplyr::select(Location,in_out, Year, N_corr) %>% 
  unique()

popsizes %>% 
  group_by(in_out) %>% 
  summarise(mean(N_corr))

mean(popsizes$N_corr)
var(popsizes$N_corr)
inout_popmod<-glmer.nb(N_corr ~ in_out + (1|Location), data = popsizes)
summary(inout_popmod)

confint.merMod(inout_popmod)

########r0 tables islands ####

isl_estimates_table<-function(tables, col_headers, parameter){
  library(tidyverse)
  
  stopifnot(length(parameter) ==1)
  
  output<-list()
  for(i in 1:length(tables)){
    
    
    table<-tables[[i]] %>% 
      dplyr::select(Location | starts_with(parameter)) %>% 
      rename(var = 2, low = 3, up = 4)
    
    table<-table %>% 
      mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>% 
      mutate(var = paste0(var, " (", low,", ", up,")")) %>% 
      dplyr::select(Location, var) %>% 
      arrange(Location) %>% 
      mutate(Location = str_to_title(Location))
    
    output[[i]] <- table
    colnames(output[[i]])[2] <- col_headers[i]
  }
  
  output<-Reduce(full_join, output)
  output<-output %>% 
    filter(Location != "Meta-Population")
}


  
#this table needs the dataframes produced in the figure above
tables <- list()

tables[[1]] <- r0_backtransformed_zNR
tables[[2]] <- zNR_single_backtransformed
tables[[3]] <- post_nometa_n
varextr <- function(x){var(x[,1])}
varextr(zNR_single_backtransformed)


col_headers <- c("Logistic, meta-population model","Logistic, single population models", "Gompertz, meta-population model")

parameter<-"r0"

(results_frame<-isl_estimates_table(tables, col_headers, parameter))

table_capt<-"Medians for posterior distributions with 90\\% quantile intervals for island-population-specific intrinsic 
growth rates $r_0$ (log individual demographic contribution at 0 adult density). Estimates come from three separate modelling 
approaches using adult density (logistic) or log adult density (Gompertz). Island-population estimates are modelled either as 
deviations from a meta-population mean $r_0$ (random effects in full meta-population models) or as fixed effects in separate models 
(single-population models)."

table_name <- "r0 table all modelling approaches"

tbl_label <- "r0_table"

library(tables)
library(knitr)
library(kableExtra)
kableExtra::kable(results_frame,booktabs = TRUE, caption = paste(table_capt),format = "latex",
                  escape = FALSE, label = paste(tbl_label),digits = 2) %>% 
  kable_paper("striped",full_width = F) %>% 
  kable_styling(font_size = 8, position = "center") %>% 
  kable_styling(latex_options="scale_down") %>% 
  save_kable(paste0("Tables empirical data/Tables output/",table_name,".tex"))

########Backcalculate variance in intercepts.####
varcors_Nm<-r2s_zNR_new %>% 
  spread_draws(sd_Location__Intercept, sd_Location__zNR_new, cor_Location__Intercept__zNR_new) %>%
  mutate(cov_isl = cor_Location__Intercept__zNR_new * sd_Location__Intercept * sd_Location__zNR_new) %>% 
  mutate(Sigma2_Isl_rmu = sd_Location__Intercept^2) %>% 
  mutate(Sigma2_Isl_gamma = sd_Location__zNR_new^2)

Nc<-data_observed %>% 
  dplyr::select(Location, Year, N_corr) %>% 
  unique() %>% 
  summarise(Nc = mean(N_corr)/40) %>% 
  as.numeric()

varcors_Nm %>% 
  mutate(Sigma2_Isl_r0 = Sigma2_Isl_rmu + Nc^2*Sigma2_Isl_gamma - 2*Nc *cov_isl) %>% 
  median_qi(Sigma2_Isl_r0, .width = 0.9)

#single pops

var(zNR_single_backtransformed$r0)

####### Gamma table all modelling approaches ####

r0_backtransformed_zNR<-r0_backtransformed_zNR %>% 
  rename(g = gamma) %>% 
  rename(g.lower = gamma.lower) %>% 
  rename(g.upper = gamma.upper)

zNR_single_backtransformed<-zNR_single_backtransformed %>% 
  rename(g = b_gamma) %>% 
  rename(g.lower = b_gamma.lower) %>% 
  rename(g.upper = b_gamma.upper)

tables <- list()

tables[[1]] <- r0_backtransformed_zNR
tables[[2]] <- zNR_single_backtransformed
tables[[3]] <- post_nometa_n
varextr <- function(x){var(x[,4])}
varextr(zNR_single_backtransformed)

col_headers <- c("Logistic, meta-population model","Logistic, single population models", "Gompertz, meta-population model")

parameter<-"g"


(results_frame<-isl_estimates_table(tables, col_headers, parameter))



table_capt<-"Medians for posterior distributions with 90\\% quantile intervals for island-population-specific strength of 
density regulation $\\gamma$ (effects of (log) adult density on individual demographic contribution (recruit production + adult survival,
negative binomial)). Estimates come from three separate modelling approaches using 
adult density (logistic) or log adult density (Gompertz). Island-population estimates are modelled either as deviations from a meta-population
mean $\\gamma$ (random effects in full meta-population models), or as fixed effects in separate models (single-population models)."

table_name <- "gamma table all modelling approaches"

tbl_label <- "gamma_table"

library(tables)
library(knitr)
library(kableExtra)
kableExtra::kable(results_frame,booktabs = TRUE, caption = paste(table_capt),format = "latex",
                  escape = FALSE, label = paste(tbl_label),digits = 2) %>% 
  kable_paper("striped",full_width = F) %>% 
  kable_styling(font_size = 8, position = "center") %>% 
  kable_styling(latex_options="scale_down") %>% 
  save_kable(paste0("Tables empirical data/Tables output/",table_name,".tex"))


#####Gamma table meta logistic with mean pop size ####

tables <- list()

tables[[1]] <- r0_backtransformed_zNR

col_headers<-"Logistic, meta-population model"

parameter <- "g"

(results_frame<-isl_estimates_table(tables, col_headers, parameter))

mean_pops<-island_N_means %>% 
  rename(`Mean density` = Mean_N_corr_isl) %>% 
  mutate(Location = str_to_title(Location))

results_frame<- results_frame %>% 
  left_join(mean_pops) %>% 
  mutate(`Mean density` = round(`Mean density`)) %>% 
  rename(Island = Location) %>% 
  arrange(`Mean density`) %>% 
  mutate(System = ifelse(in_out == 0, "Outer", "Inner")) %>% 
  relocate(System, .after = Island) %>% 
  dplyr::select(-in_out)
           


table_capt<-"Medians for posterior distributions with 90\\% quantile intervals for island-population-specific strength of 
density regulation $\\gamma$ (effects of adult density on individual demographic contribution (recruit production + adult survival,
negative binomial)) and mean island-population density."

table_name <- "gamma density table"

tbl_label <- "gamma_mean_density_table"

library(tables)
library(knitr)
library(kableExtra)
kableExtra::kable(results_frame,booktabs = TRUE, caption = paste(table_capt),format = "latex",
                  escape = FALSE, label = paste(tbl_label),digits = 2) %>% 
  kable_paper("striped",full_width = F) %>% 
  kable_styling(font_size = 8, position = "center") %>% 
  kable_styling(latex_options="scale_down") %>% 
  save_kable(paste0("Tables empirical data/Tables output/",table_name,".tex"))


#######Fitness at the mean density tables islands ####

post_nb_full_Nm_nometa<-post_nb_full_Nm %>% 
  filter(Island != "metapopulation")


n_params<-r2s_n_new %>% 
  spread_draws(b_Intercept, b_n_new, b_in_out, `b_n_new:in_out`, r_Location[Location, variable])

n_params<-n_params %>% 
  pivot_wider(names_from = variable, values_from = r_Location)

n_params$Location <- gsub("\\.", " ", n_params$Location)
unique(n_params$Location)

n_params <-n_params %>% 
  left_join(island_N_means)

n_params<-n_params %>% 
  mutate(r0 = ifelse(in_out == 0,b_Intercept + Intercept -log(2), b_Intercept +b_in_out + Intercept -log(2))) %>% 
  mutate(g = ifelse(in_out == 0, b_n_new +  n_new, b_n_new + `b_n_new:in_out` + n_new))

n_params %>% 
  group_by(Location) %>% 
  median_qi(Intercept)
n_params %>% 
  ungroup() %>% 
  median_qi(b_Intercept, b_in_out)

0.983 + -0.0757 + 0.0216-log(2)
post_nb_full
n_params %>% 
  group_by(Location) %>%
  median_qi(r0,g)
0.254 + -0.165 *log(148)
n_rmu_tab<-n_params %>% 
  group_by(Location) %>% 
  mutate(r0 = r0 +  g*log(Mean_N_corr_isl)) %>% 
  median_qi(r0,g, .width=0.9) %>% 
  ungroup()

post_nb_full_Nm_nometa<-post_nb_full_Nm_nometa %>% 
  rename(Location = Island)
post_frame_single<-post_frame_single %>% 
  rename(Location = Island)
tables <- list()

tables[[1]] <- post_nb_full_Nm_nometa
tables[[2]] <- post_frame_single
tables[[3]] <- n_rmu_tab
varextr <- function(x){var(x[,5])}

lapply(tables, varextr)

col_headers <- c("Logistic, meta-population model","Logistic, single population models", "Gompertz, meta-population model")

parameter<-"r0"

(results_frame<-isl_estimates_table(tables, col_headers, parameter))

table_capt<-"Medians for posterior distributions with 90\\% quantile intervals for island-population-specific 
log individual demographic contribution at island-population-specific mean adult density $r_{\\mu}$. 
Estimates come from three separate modelling approaches using adult density (logistic) or log adult density (Gompertz). 
Island-population estimates are modelled either as deviations from a meta-population
mean $r_{\\mu}$ (random effects in full meta-population models) or as fixed effects in separate models (single-population models)."

table_name <- "mean density fitness table all modelling approaches"

tbl_label <- "mean_density_fitness_table"

library(tables)
library(knitr)
library(kableExtra)
kableExtra::kable(results_frame,booktabs = TRUE, caption = paste(table_capt),format = "latex",
                  escape = FALSE, label = paste(tbl_label),digits = 2) %>% 
  kable_paper("striped",full_width = F) %>% 
  kable_styling(font_size = 8, position = "center") %>% 
  kable_styling(latex_options="scale_down") %>% 
  save_kable(paste0("Tables empirical data/Tables output/",table_name,".tex"))


###### Vital rates used for elasticities table ####

nest_to_rec_IPM_genquant_block<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")

isl_ids<-data_observed %>% 
  dplyr::select(Location, Island, in_out) %>% 
  filter(Island != 33) %>% 
  unique() %>% 
  mutate(Island_stan = as.numeric(as.factor(Island)))

nest_to_rec_IPM_genquant_block %>% 
  spread_draws(f_simple[Island_stan]) %>% 
  ungroup() %>% 
  median_qi(f_simple, .width=0.9)

f_simple<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(f_simple[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(f_simple, .width=0.9) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, f_simple, .lower, .upper)


f_simple<-f_simple %>% 
  left_join(isl_ids) %>% 
  mutate("$F$" = paste(round(f_simple,2), " (", round(.lower,2), ",", round(.upper,2),")", sep = "")) %>% 
  dplyr::select(Location,in_out, "$F$") %>% 
  arrange(Location)

f_simple_in_out<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(f_simple[Island_stan]) 

f_simple_in_out<-f_simple_in_out %>% 
  left_join(isl_ids) %>% 
  group_by(in_out) %>% 
  median_qi(f_simple,.width = 0.9) %>% 
  mutate("$F$" = paste(round(f_simple,2), " (", round(.lower,2), ",", round(.upper,2),")", sep = "")) %>% 
  mutate(Location = NA) %>% 
  dplyr::select(Location,in_out, "$F$")

f_simple<-bind_rows(f_simple_in_out, f_simple)

s_n<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(s_n[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(s_n, .width=0.9) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, s_n, .lower, .upper)

s_n<-s_n %>% 
  left_join(isl_ids) %>% 
  mutate("$\\phi_{nr}$" = paste(round(s_n,2), " (", round(.lower,2), ",", round(.upper,2),")", sep = "")) %>% 
  dplyr::select(Location,in_out, "$\\phi_{nr}$")  %>% 
  arrange(Location)

s_a<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(s_a[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(s_a, .width=0.9) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, s_a, .lower, .upper)

s_a<-s_a %>% 
  left_join(isl_ids) %>% 
  mutate("$\\phi_{ad}$" = paste(round(s_a,2), " (", round(.lower,2), ",", round(.upper,2),")", sep = "")) %>% 
  dplyr::select(Location,in_out, "$\\phi_{ad}$")  %>% 
  arrange(Location)



outer_ad<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(mu_phi)%>%
  mutate(outer_ad = invlogit(mu_phi)) %>% 
  median_qi(outer_ad, .width =0.9)

outer_ad<-outer_ad %>% 
  mutate("$\\phi_{ad}$" = paste(round(outer_ad,2), " (",round(.lower,2), ", ", round(.upper,2), ")", sep ="")) %>% 
  mutate(in_out = 0) %>% 
  dplyr::select(in_out, "$\\phi_{ad}$")

inner_ad<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(mu_phi, B_in_out)%>%
  mutate(inner_ad = invlogit(mu_phi + B_in_out)) %>% 
  median_qi(inner_ad, .width =0.9)

inner_ad<-inner_ad %>% 
  mutate("$\\phi_{ad}$" = paste(round(inner_ad,2), " (",round(.lower,2), ", ", round(.upper,2), ")", sep ="")) %>% 
  mutate(in_out = 1) %>% 
  dplyr::select(in_out, "$\\phi_{ad}$")


ad_inout<-bind_rows(outer_ad, inner_ad)

ad_inout<-ad_inout %>% 
  mutate(Location =NA) %>% 
  relocate(Location, .before = in_out)

outer_nr<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(mu_phi,mu_phi_rec) %>% 
  mutate(outer_nr = invlogit(mu_phi +mu_phi_rec)) %>% 
  median_qi(outer_nr, .width =0.9) %>% 
  dplyr::select(outer_nr, .lower, .upper)

outer_nr<-outer_nr %>% 
  mutate("$\\phi_{nr}$" = paste(round(outer_nr,2), " (",round(.lower,2), ", ", round(.upper,2), ")", sep ="")) %>% 
  mutate(in_out = 0) %>% 
  dplyr::select(in_out, "$\\phi_{nr}$")

inner_nr <-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(mu_phi,mu_phi_rec, B_in_out_rec) %>% 
  mutate(inner_nr = invlogit(mu_phi +mu_phi_rec + B_in_out_rec)) %>% 
  median_qi(inner_nr, .width =0.9)

inner_nr<-inner_nr %>% 
  mutate("$\\phi_{nr}$" = paste(round(inner_nr,2), " (",round(.lower,2), ", ", round(.upper,2), ")", sep ="")) %>% 
  mutate(in_out = 1) %>% 
  dplyr::select(in_out, "$\\phi_{nr}$")

nr_inout<-bind_rows(outer_nr, inner_nr)

nr_inout<-nr_inout %>% 
  mutate(Location =NA)%>% 
  relocate(Location, .before = in_out)

s_a <- bind_rows(ad_inout, s_a)
s_n <-bind_rows(nr_inout, s_n)

proj_matrix_vital_rates <- f_simple %>% 
  left_join(s_n) %>% 
  left_join(s_a) %>% 
  mutate(Location = str_to_title(Location)) %>% 
  mutate(in_out = ifelse(in_out == 0, "Outer non-farm", "Inner farm")) %>% 
  rename(`Habitat system` = in_out)


# To make it simpler we create captions in footnotes in Overleaf
library(xtable)

output_filename <- "proj_matrix_vital_rates.tex"

xtable_result <- xtable(proj_matrix_vital_rates, 
                        floating = FALSE, 
                        escape = FALSE, 
                        booktabs = TRUE)

print(
  xtable_result, 
  file = paste0("Tables empirical data/Tables output/", output_filename), 
  include.rownames = FALSE, 
  floating = FALSE, #setting floating to false removes \begin{table} and \end{table} from the file
  sanitize.text.function = identity, 
  sanitize.colnames.function = function(x) { paste0("\\multicolumn{1}{l}{", x, "}") },
)

######Elasticity table####

nest_to_rec_IPM_genquant_block<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")

inout<-data_observed %>% 
  dplyr::select(Location,Island, in_out) %>% 
  filter(Island != 33) %>%
  mutate(Island_stan = as.numeric(as.factor(Island))) %>% 
  dplyr::select(-Island) %>% 
  unique() %>% 
  arrange(Island_stan)

elast_in_out<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(elasticity[Island_stan, row, col])

elast_in_out<-elast_in_out %>% 
  left_join(inout)


elast_in_out<-elast_in_out %>% 
  group_by(in_out, row, col) %>% 
  median_qi(elasticity=round(elasticity,2), .width = 0.9) %>% 
  ungroup()

elast_in_out$elasticity <- paste(elast_in_out$elasticity, " (", elast_in_out$.lower, ", ", elast_in_out$.upper, ")", sep ="")

elast_in_out<-elast_in_out %>%
  mutate(Location = NA) %>% 
  relocate(Location, .before = in_out)

elast_frame<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(elasticity[Island_stan, row, col]) %>% 
  group_by(Island_stan ,row, col) %>% 
  median_qi(elasticity=round(elasticity,2), .width = 0.9) %>% 
  ungroup() 

elast_frame$elasticity <- paste(elast_frame$elasticity, " (", elast_frame$.lower, ", ", elast_frame$.upper, ")", sep ="")

elast_frame<-elast_frame %>% 
  left_join(inout) %>% 
  arrange(Location)

elast_frame<-elast_frame %>% 
  dplyr::select( Location, in_out,row,col, elasticity,.lower, .upper, .width, .point, .interval)



elast_frame<-bind_rows(elast_in_out, elast_frame)

elast_frame<-elast_frame %>% 
  unite("rowcol", row,col, sep ="_") %>% 
  dplyr::select(Location, in_out, elasticity, rowcol) %>% 
  pivot_wider(values_from = elasticity, names_from = rowcol) %>% 
  rename("$F$" = `1_2`) %>% 
  rename("$\\phi_{nr}$" =`2_1`) %>% 
  rename( "$\\phi_{ad}$" = `2_2`) %>% 
  rename(`Habitat system`= in_out) %>% 
  dplyr::select(-`1_1`) %>% 
  mutate(Location = str_to_title(Location))

elast_frame$`Habitat system`[elast_frame$`Habitat system` == 0] <- "Outer non-farm"
elast_frame$`Habitat system`[elast_frame$`Habitat system` == 1] <- "Inner farm"

elast_frame



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



#####Cumulative estimated nestlings, recruits and adults and cumulative observed recruits and adults table ####

nest_to_rec_IPM_genquant_block<-readRDS("Workspace backup/nest_to_rec_IPM_genquant_block.rds")
demo4<-readRDS("Workspace backup/demo4.rds")

##Number of nestlings per habitat system
nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Nest_corr_F[Island_stan]) %>% 
  left_join(isl_ids) %>% 
  group_by(.draw, in_out) %>% 
  summarise(Nest_corr_F_sum = sum(Nest_corr_F), .groups = "drop") %>% 
  group_by(in_out) %>% 
  median_qi(Nest_corr_F_sum, .width = 0.9) %>% 
  ungroup()

Nest_corr_f<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Nest_corr_F[Island_stan]) %>% 
  group_by(Island_stan) %>% 
  median_qi(Nest_corr_F, .width = 0.9) %>% 
  ungroup() %>% 
  dplyr::select(Island_stan, Nest_corr_F, .lower, .upper, .width)

Nest_corr_f<-Nest_corr_f %>% 
  left_join(isl_ids)

Nest_corr_f %>% 
  group_by(in_out) %>% 
  summarise(sum(Nest_corr_F))

Nest_corr_summary<-Nest_corr_f %>% 
  mutate(`Est. nest.` = paste(round(Nest_corr_F), " (", round(.lower),",", round(.upper), ")", sep = "")) %>% 
  dplyr::select(Location, in_out, `Est. nest.`) %>% 
  arrange(Location)

##Total recruits per inner and outer
nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Rec_corr[Island_stan, Year_stan]) %>% 
  left_join(isl_ids) %>% 
  group_by(.draw, in_out) %>% 
  summarise(Rec_corr_sum = sum(Rec_corr), .groups = "drop") %>% 
  group_by(in_out) %>% 
  median_qi(Rec_corr_sum, .width = 0.9) %>% 
  ungroup()

Rec_corr_summary<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(Rec_corr[Island_stan, Year_stan]) %>% 
  group_by(.draw, Island_stan) %>% 
  summarise(Rec_corr_sum = sum(Rec_corr), .groups = "drop") %>% 
  group_by(Island_stan) %>% 
  median_qi(Rec_corr_sum, .width = 0.9) %>% 
  ungroup()


Rec_corr_summary<-Rec_corr_summary %>% 
  left_join(isl_ids) %>% 
  mutate(`Est. rec.` = paste(round(Rec_corr_sum), " (", round(.lower),",", round(.upper), ")", sep = "")) %>% 
  dplyr::select(Location, in_out, `Est. rec.`) %>% 
  arrange(Location)

Obs_rec_summary<-demo4 %>% 
  group_by(Location) %>% 
  summarise( `Obs. rec.` =sum(Rec_obs))

Obs_ad_summary<-demo4 %>% 
  group_by(Location) %>% 
  summarise( `Obs. ad.` =sum(N_obs))

#adults per habitat system
nest_to_rec_IPM_genquant_block %>% 
  spread_draws(N_corr[Island_stan, Year_stan]) %>% 
  left_join(isl_ids) %>% 
  group_by(.draw, in_out) %>% 
  summarise(N_corr_sum = sum(N_corr), .groups = "drop") %>% 
  group_by(in_out) %>% 
  median_qi(N_corr_sum, .width = 0.9) %>% 
  ungroup()


N_corr_summary<-nest_to_rec_IPM_genquant_block %>% 
  spread_draws(N_corr[Island_stan, Year_stan]) %>% 
  group_by(.draw, Island_stan) %>% 
  summarise(N_corr_sum = sum(N_corr), .groups = "drop") %>% 
  group_by(Island_stan) %>% 
  median_qi(N_corr_sum, .width = 0.9) %>% 
  ungroup()


N_corr_summary<-N_corr_summary %>% 
  left_join(isl_ids) %>% 
  mutate(`Est. ad.` = paste(round(N_corr_sum), " (", round(.lower),",", round(.upper), ")", sep = "")) %>% 
  dplyr::select(Location, in_out, `Est. ad.`) %>% 
  arrange(Location)


individuals_obs_est<-Nest_corr_summary %>% 
  left_join(Obs_rec_summary) %>% 
  left_join(Rec_corr_summary) %>% 
  left_join(Obs_ad_summary) %>% 
  left_join(N_corr_summary) %>% 
  mutate(Location = str_to_title(Location)) %>% 
  mutate(`Habitat system` = ifelse(in_out == 0, "Outer non-farm system", "Inner farm system")) %>% 
  dplyr::select(Location, `Habitat system`, `Est. nest.`, `Obs. rec.`, `Est. rec.`, `Obs. ad.`, `Est. ad.`)


output_filename <- "cumulative_nestling_recruits_adults_table.tex"

xtable_result <- xtable(individuals_obs_est, 
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


