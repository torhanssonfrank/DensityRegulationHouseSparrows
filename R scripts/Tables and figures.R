
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
data_observed<-read.table("Data/cleaned_density_dependence_data.txt", sep=";", header=TRUE)

##Load latex variable names ####

paramlatex<-read_xlsx(path="Tables empirical data/Parameter to latex.xlsx")


##Negative binomial log N table, no individual slopes####
bernt_n_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))
n_table<-brm2table_short(bernt_n_nb_brm_noindslope, type ="median" , width =0.9)

table_capt<-"Medians for posterior distributions for the effects of log density 
(log adult density) and other sources of variation 
on individual demographic contribution (recruit production + adult survival, with
90\\% credible intervals in parentheses."

tables <- list()

tables[[1]] <- n_table

col_headers<-c("\\textbf{Recruit prod. + Survival}")
tbl_label <- "n_nb_tbl"
table_name <- "results_fixed_n_noindslope"
table_footnote <-NULL
sigfigs = 2
n_fixef = 8
fit_table<-table_wrangler_smart(tables,col_headers,sigfigs,n_fixef, paramlatex)

res_tables(fit_table, table_capt, tbl_label,table_footnote, table_name, n_fixef)

##r2S, recruits and adult survival CMR N_m table####
bernt_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))
rec_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/rec_zNR_nb_brm_noindslope.rds"))
CMR_adult_survival<-readRDS("Workspace backup/CMR_adults_interacs_old_Paul_age_sq_zNR.rds")
mean_age<-mean(bernt_zNR_nb_brm_noindslope$data$Least_age)


width =0.9
bernt_zNR_nb_brm_noindslope_table <-brm2table_short(bernt_zNR_nb_brm_noindslope, type ="median", width)
rec_zNR_nb_brm_noindslope_table<-brm2table_short(rec_zNR_nb_brm_noindslope, type = "median", width)

reverse_inner_outer="no"

CMR_adult_table<-CMR2table_adults(CMR_adult_survival,paramlatex,reverse_inner_outer, width)



roundinvlogit<-function(x){round(invlogit(x), 2)
}
mu_p<-CMR_adult_table %>% 
  filter(Parameter == "mu_p") %>% 
  mutate(across(where(is.numeric), roundinvlogit))

Sigma2_YI_p<-CMR_adult_table %>% 
  filter(Parameter == "Sigma2_YI_p") %>% 
  mutate(across(where(is.numeric), ~round(.x,2)))


table_capt<-"Medians for posterior distributions for the effects of density 
(mean-centred and scaled adult density) and other sources of variation 
on individual demographic contribution (recruit production + adult survival, 
negative binomial) and its components: recruit production (negative binomial) 
and adult survival (binomial, mark-recapture) with
90\\% credible intervals in parentheses."

tbl_label <- "Nm_nb_tbl"
table_name <- "results_fixed_nb_brm_noindslope"
table_footnote <- paste0("The average recapture probability was p = ",mu_p$Estimate, " (90\\% CI = ", mu_p$.lower, ", ", mu_p$.upper, ") with a variance \\\ among years within populations on the latent scale of ", Sigma2_YI_p$Estimate, " (90\\% CI = ", Sigma2_YI_p$.lower, ", ", Sigma2_YI_p$.upper, "). Island, Island:density 
and Island cov represent the (co)variances of island-specific intercepts and slopes 
for the effect of density. Individual and Island-year represent variances in individual and Island-year 
specific intercepts. The unexplained variance in adult survival is not presented since 
it is captured by the Bernoulli process.")


tables <- list()

tables[[1]] <- bernt_zNR_nb_brm_noindslope_table
tables[[2]] <- rec_zNR_nb_brm_noindslope_table
tables[[3]] <- CMR_adult_table

col_headers<-c("\\textbf{Demographic contribution}", "\\textbf{Recruits}", "\\textbf{Adult survival}")
sigfigs = 2
n_fixef = 8
fit_table<-table_wrangler_smart(tables,col_headers,sigfigs,n_fixef, paramlatex)
fit_table<-fit_table %>% 
  filter(str_detect( `\\textbf{Parameter}`,"Mean recapture probability|Sigma2_YI_p", negate = T))
fit_table

#results_frame = fit_table
footn <- list(pos=list(0), command= NULL)
footn$pos[[1]] <- c(nrow(fit_table))
footn$command <- c(paste0("\\hline\\multicolumn{4}{l}{Footnote: See Appendix S6 for recapture probabilities. Island, Island:density and Island cov}\\\\
                          \\multicolumn{4}{l}{represent the (co)variances of island-specific intercepts and slopes for the effect of density.}\\\\
                          \\multicolumn{4}{l}{Individual and Island-year represent variances in individual and Island-year specific intercepts.}\\\\
                          \\multicolumn{4}{l}{The unexplained variance in adult survival is not presented since it is captured by the Bernoulli}\\\\
                          \\multicolumn{4}{l}{process.}", "\n", sep = ""))


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

#res_tables(results_frame, table_capt, tbl_label,table_footnote, table_name, n_fixef)

## Juvenile survival CMR and nestling production ####

CMR_juv<-readRDS("Workspace backup/CMR_juvenile_survival_Tor_zNR.rds")
nestlings_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/nestlings_zNR_nb_brm_noindslope.rds")) 

width = 0.9
nestling_prod_res <- brm2table_short(nestlings_zNR_nb_brm_noindslope, type = "median",width)

reverse_inner_outer="yes"
add_difference_to_intercept="yes"

juv_res<-CMR2table_juveniles(CMR_juv,paramlatex,reverse_inner_outer,add_difference_to_intercept, width)


(nest2fledge_res<-juv_res[[1]])
(fledge2rec_res<-juv_res[[2]])
nestling_prod_res
#clutch_res

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
footn$command <- c(paste0(
  "\\hline","\\multicolumn{4}{l}{Footnote: Nestling and fledgling survival was modelled in a mark-recapture framework. The}", "\\\\", 
  "\\multicolumn{4}{l}{average fledgling recapture was p = ", mu_p_fledge$Estimate,  " (90\\% CI = ", mu_p_fledge$.lower, ", ", mu_p_fledge$.upper, ") with a variance among years}", "\\\\", 
  "\\multicolumn{4}{l}{within populations on the latent scale of ", Sigma2_YI_p_fledge$Estimate, " (90\\% CI = ", Sigma2_YI_p_fledge$.lower, ", ", Sigma2_YI_p_fledge$.upper, "). The average recruit recapture}","\\\\",
  "\\multicolumn{4}{l}{probability was p = ",mu_rec_p$Estimate, " (90\\% CI = ", mu_rec_p$.lower, ", ", mu_rec_p$.upper, ") with a variance among years within populations}", "\\\\",
  "\\multicolumn{4}{l}{on the latent scale of ", Sigma2_YI_p_rec$Estimate," (90\\% CI = ", Sigma2_YI_p_rec$.lower, ", ", Sigma2_YI_p_rec$.upper, "). Island, Island:density and Island cov represent}", "\\\\", 
  "\\multicolumn{4}{l}{the (co)variances of island-specific intercepts and slopes for the effect of density. Individual and}", "\\\\", 
  "\\multicolumn{4}{l}{Island-year represent variances in individual and Island-year specific intercepts.The unexplained}", "\\\\",
  "\\multicolumn{4}{l}{variance in Bernoulli models is not presented since it is captured by the Bernoulli process.}"))

#\\hline,\\multicolumn{4}{l}{Footnote: See Appendix S7 for recapture probabilities. Island, Island:density and Island cov}\\\multicolumn{4}{l}{represent the (co)variances of island-specific intercepts and slopes for the effect of density. Individual}\\\multicolumn{4}{l}{and Island-year represent variances in individual and Island-year specific intercepts. The unexplained }\\\multicolumn{4}{l}{variance in Bernoulli models is not presented since it is captured by the Bernoulli process.}


#add.to.row = footn, hline.after = c(-1,0),
print(xtable(fit_table,type = "latex", caption = table_capt, label = tbl_label,floating=FALSE, escape = F,booktabs = T), file = paste0("Tables empirical data/Tables output/",table_name,".tex"),caption.placement = "top", include.rownames=FALSE,
      add.to.row = footn, hline.after = c(-1,0),
      sanitize.text.function = identity, sanitize.colnames.function=function(x){paste0("\\multicolumn{1}{l}{",x,"}")},size="\\fontsize{12pt}{20pt}\\selectfont") #multicolumn{1}{l} the l left-aligns the col headers

#res_tables(fit_table, table_capt, tbl_label, table_footnote, table_name,n_fixef)
## r2s Negative binomial log N table, and r2s mean centred table####
bernt_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))
zNR_table <-brm2table_short(bernt_zNR_nb_brm_noindslope, type ="median")
bernt_n_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))
logn_table<-brm2table_short(bernt_n_nb_brm_noindslope, type ="median")
VarCorr(bernt_n_nb_brm_noindslope)$ID

table_capt<-"Results of the random regression models for density regulation modelled with mean-centred population size and log population size on recruit production and survival."
tbl_label <- "n_mcN_nb_tbl"
table_name <- "results_fixed_nb_n_mcN_noindslope"
tables <- list()

tables[[1]] <- zNR_table
tables[[2]] <- logn_table
col_headers<-c("$N_{z{\\mu}}$", "$n$")
fit_table<-table_wrangler_smart(tables,col_headers, paramlatex)
fit_table

n_fixef = 7

res_tables(fit_table, table_capt, tbl_label, table_name, n_fixef)


# Island predicted slopes at experienced pop sizes  ####
##log N ####
#r2s

density_variable = "n"
fitness = "r2s"
divide_intercepts= "yes"
average_metapars ="yes"

bernt_n_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))
brm2table_short(bernt_n_nb_brm_noindslope, type ="median" , width =0.9)

post_nb_full<-brm_blups_islands(bernt_n_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)


post_frame= post_nb_full

popsize_frame <- data_observed %>% 
  dplyr::select(Location, n) %>% 
  rename(popsize = n) %>% 
  rename(Island = Location)
mean_centred = "no"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5


island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)
##Save v on n and w on N island plots to jpeg####

distribution = "nb"
jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100)
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5)
# 2. Create a plot

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()


## Island slopes mean centred ####
#zmc_N (zNR)

density_variable = "zmc_N"
fitness = "r2s"
divide_intercepts = "yes"
average_metapars ="yes"
mean_centred = "yes"
bernt_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))
summary(bernt_zNR_nb_brm_noindslope)
brm2table_short(bernt_zNR_nb_brm_noindslope, type="median", width = 0.9)

post_nb_full<-brm_blups_islands(bernt_zNR_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)
post_nb_full

source("Functions/brm2table_short.r")

brm2table_short(bernt_zNR_nb_brm_noindslope, type = "median", width = 0.9)

post_frame= post_nb_full

popsize_frame <- data_observed %>% 
  dplyr::select(Location, zNR) %>% 
  rename(popsize = zNR) %>% 
  rename(Island = Location)



distribution = "nb"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

#Save v on n and w on N island plots to jpeg
jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100) 
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_islands_meta_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5)
# 2. Create a plot

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()


## Islands modelled as separate populations ####

#r2s zNR

bernt_zNR_single_island_mods<-readRDS(file("Workspace backup/bernt_zNR_single_island_mods.rds"))
bernt_zNR_single_island_mods[[1]]
divide_intercepts = "yes"

post_frame= get_single_island_estimates(bernt_zNR_single_island_mods,divide_intercepts, summary = TRUE)


popsize_frame <- data_observed %>% 
  dplyr::select(Location, zNR) %>% 
  rename(popsize = zNR) %>% 
  rename(Island = Location)


mean_centred ="yes"

density_variable ="zmc_N"

prob = 1
v_ymin = -1.2
v_ymax =0.35
w_ymin = 0.35
w_ymax = 1.5
fitness = "r2s"
island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

distribution = "nb"

jpeg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_single_islands_v_",density_variable,"_wN_predicted.jpeg", sep = ""),
     width = 2000, height = 1100) 
#svg(paste("Figures empirical data/",distribution,"_prob_",prob,"_",fitness, "_", "_single_islands_v_",density_variable,"_wN_predicted.svg", sep = ""),
#    width = 30, height = 16.5) 
# 2. Create a plot

island_slopes(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness)

# Close the jpeg file
dev.off()


#Gamma table, single islands, metapop zNR, metapop n ####
bernt_zNR_single_island_mods<-readRDS(file("Workspace backup/bernt_zNR_single_island_mods.rds"))
divide_intercepts = "yes"

post_single_zNR= get_single_island_estimates(bernt_zNR_single_island_mods,divide_intercepts, summary = TRUE)




density_variable = "zmc_N"
fitness = "r2s"
average_metapars ="yes"
mean_centred = "yes"
bernt_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))



post_metapop_zNR<-brm_blups_islands(bernt_zNR_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)



density_variable = "n"
mean_centred = "no"

bernt_n_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))


post_metapop_n<-brm_blups_islands(bernt_n_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)

tables <- list()

tables[[1]] <- post_metapop_zNR
tables[[2]] <- post_single_zNR
tables[[3]] <- post_metapop_n
varextr <- function(x){var(x[,5])}
varextr(post_single_zNR)
lapply(tables, varextr)

col_headers <- c("Logistic, meta-population model","Logistic, single population models", "Gompertz, meta-population model")

parameter<-"g"

isl_estimates_table<-function(tables, col_headers, parameter){
  library(tidyverse)
  
  stopifnot(length(parameter) ==1)
  
  output<-list()
  for(i in 1:length(tables)){
    
    
    table<-tables[[i]] %>% 
      dplyr::select(Island | starts_with(parameter)) %>% 
      rename(var = 2, low = 3, up = 4)
    
    table<-table %>% 
      mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>% 
      mutate(var = paste0(var, " (", low,", ", up,")")) %>% 
      dplyr::select(Island, var) %>% 
      arrange(Island) %>% 
      mutate(Island = str_to_title(Island))
    
    output[[i]] <- table
    colnames(output[[i]])[2] <- col_headers[i]
  }
  
  output<-Reduce(full_join, output)
  output<-output %>% 
    filter(Island != "Meta-Population")
}

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


tables <- list()



tables[[1]] <- post_metapop_zNR

col_headers<-"Logistic, meta-population model"

parameter <- "g"

(results_frame<-isl_estimates_table(tables, col_headers, parameter))

mean_pops<-island_N_means %>% 
  rename(Island = Location) %>% 
  rename(`Mean density` = Mean_N_isl) %>% 
  mutate(Island = str_to_title(Island))

results_frame<- results_frame %>% 
  left_join(mean_pops) %>% 
  mutate(`Mean density` = round(`Mean density`)) %>% 
  arrange(`Mean density`)


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

##Life stage flow chart figure ####

nestlings_zNR_nb_brm_noindslope %>% 
  get_variables() %>% 
  head(n=12)

mean_age<-mean(nestlings_zNR_nb_brm_noindslope$data$Least_age)

nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_zNR, `b_zNR:in_out`) %>% 
  mutate(meta_g = b_zNR + `b_zNR:in_out`*0.5) %>% 
  median_qi(meta_g, .width=0.9)


nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_zNR, `b_zNR:in_out`) %>% 
  mutate(inner_g = b_zNR + `b_zNR:in_out`) %>% 
  median_qi(round(inner_g,2), .width=0.9)

nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_zNR) %>% 
  mutate(outer_g = b_zNR) %>% 
  median_qi(round(outer_g,2), .width=0.9)

nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_Intercept, b_in_out, b_Least_age, b_Least_age_sq) %>% 
  mutate(meta_prod = b_Intercept + b_in_out*0.5+ b_Least_age*mean_age 
         + b_Least_age_sq*mean_age^2) %>% 
  median_qi(exp(meta_prod-log(2)), .width=0.9)


nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_Intercept, b_in_out, b_Least_age, b_Least_age_sq) %>% 
  mutate(inner_prod = b_Intercept + b_in_out+ b_Least_age*mean_age 
         + b_Least_age_sq*mean_age^2) %>% 
  median_qi(exp(inner_prod-log(2)), .width=0.9)

nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_Intercept, b_in_out, b_Least_age, b_Least_age_sq,
               b_sex, `b_Least_age:sex`) %>% 
  mutate(outer_prod = b_Intercept + b_Least_age*mean_age 
         + b_Least_age_sq*mean_age^2
  ) %>% 
  median_qi(exp(outer_prod-log(2)), .width=0.9)

nestlings_zNR_nb_brm_noindslope %>%
  spread_draws(b_Intercept, b_in_out) %>% 
  mutate(meta_prod = exp(b_Intercept + b_in_out*0.5)/2) %>% 
  median_qi(meta_prod, .width=0.9)


CMR_juv %>% 
  get_variables() %>% 
  head(n = 12)

#inner is 1 in this model.
CMR_juv %>% 
  spread_draws(mu_phi) %>% 
  median_qi(outer_fledge_intercept = round(invlogit(mu_phi),2), .width = 0.9) 
CMR_juv %>% 
  spread_draws(mu_phi, phi_fledge_interac_inout) %>% 
  median_qi(inner_fledge_intercept = round(invlogit(mu_phi + phi_fledge_interac_inout),2), .width = 0.9)


CMR_juv %>% 
  spread_draws(mu_phi, phi_fledge_interac_inout) %>% 
  mutate(inner_fledge_intercept = mu_phi + phi_fledge_interac_inout) %>% 
  mutate(prop_fledge_phi = invlogit(inner_fledge_intercept)/invlogit(mu_phi)) %>% 
  median_qi(prop_fledge_phi, .width =0.9)

CMR_juv %>% 
  spread_draws(mu_phi, phi_fledge_interac_inout) %>% 
  mutate(inner_fledge_intercept = mu_phi + phi_fledge_interac_inout) %>% 
  mutate(outer_fledge_diff=  mu_phi - inner_fledge_intercept) %>% 
  median_qi(outer_fledge_diff, .width = 0.9)

CMR_juv %>% 
  spread_draws(phi_fledge_interac_inout) %>%
  median_qi(phi_fledge_interac_inout, .width = 0.9)

CMR_juv %>% 
  spread_draws(gamma) %>% 
  median_qi(outer_fledge_gamma = gamma, .width = 0.9) %>% 
  mutate(outer_fledge_gamma = round(outer_fledge_gamma,2))



CMR_juv %>% 
  spread_draws(gamma,gamma_fledge_interac_inout ) %>% 
  mutate(inner_gamma= gamma + gamma_fledge_interac_inout) %>% 
  mutate(outer_fledge_gamma_diff = gamma-  inner_gamma) %>% 
  median_qi(outer_fledge_gamma_diff, .width = 0.9) %>% 
  mutate(outer_fledge_gamma_diff = round(outer_fledge_gamma_diff,6),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(gamma, gamma_fledge_interac_inout ) %>% 
  median_qi(inner_fledge_gamma = gamma + gamma_fledge_interac_inout , .width = 0.9) %>% 
  mutate(inner_fledge_gamma = round(inner_fledge_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))


#latent intercepts recruits

CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi) %>% 
  median_qi(outer_rec_phi =  mu_phi + mu_rec_phi, .width = 0.9) %>% 
  mutate(outer_rec_phi = round(outer_rec_phi,3),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout) %>% 
  median_qi(inner_rec_phi = mu_phi+mu_rec_phi + phi_rec_interac_inout, .width = 0.9) %>% 
  mutate(inner_rec_phi = round(inner_rec_phi,3),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout) %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi)) %>%
  mutate(prop_rec_phi = outer_rec_phi/inner_rec_phi) %>% 
  median_qi(prop_rec_phi, .width = 0.9) %>% 
  mutate(prop_rec_phi = (1-prop_rec_phi)*100) %>% 
  mutate(.lower = (1-.lower)*100) %>% 
  mutate(.upper = (1-.upper)*100)


CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout, phi_fledge_interac_inout) %>% 
  mutate(outer_rec_phi = mu_phi + mu_rec_phi) %>% 
  mutate(inner_rec_phi = mu_phi+mu_rec_phi + phi_rec_interac_inout) %>% 
  mutate(inner_mu_phi = mu_phi + phi_fledge_interac_inout) %>% 
  median_qi(inner_rec_phi_diff = inner_rec_phi - inner_mu_phi, .width = 0.9) %>% #inner_rec_phi_diff should be diff to fledglings in inner, not outer.
  mutate(inner_rec_phi_diff = round(inner_rec_phi_diff,3),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))


CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout, phi_fledge_interac_inout) %>% 
  mutate(outer_rec_phi = mu_phi + mu_rec_phi) %>% 
  mutate(inner_rec_phi = mu_phi+mu_rec_phi + phi_rec_interac_inout) %>% 
  mutate(inner_mu_phi = mu_phi + phi_fledge_interac_inout) %>% 
  mutate(inner_rec_phi_diff_to_inner_fledge = inner_rec_phi - inner_mu_phi) %>% 
  median_qi(inner_rec_phi_diff_to_inner_fledge, .width =0.9) %>% 
  mutate(inner_rec_phi_diff_to_inner_fledge = round(inner_rec_phi_diff_to_inner_fledge,3),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))


CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout) %>% 
  mutate(outer_rec_phi = mu_phi + mu_rec_phi) %>% 
  mutate(inner_rec_phi = mu_phi+mu_rec_phi + phi_rec_interac_inout) %>% 
  median_qi(inner_rec_phi_diff = inner_rec_phi - outer_rec_phi, .width = 0.9) %>% 
  mutate(inner_rec_phi_diff = round(inner_rec_phi_diff,3),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

juv_res
#prob intercepts recruits
CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi) %>% 
  median_qi(outer_rec_phi =  invlogit(mu_phi + mu_rec_phi), .width = 0.9) %>% 
  mutate(outer_rec_phi = round(outer_rec_phi,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(phi_rec_interac_inout) %>% 
  median_qi(inner_rec_phi_diff = phi_rec_interac_inout, .width = 0.9) %>% 
  mutate(inner_rec_phi_diff = round(inner_rec_phi_diff,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(mu_phi, mu_rec_phi,phi_rec_interac_inout) %>% 
  median_qi(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout), .width = 0.9) %>% 
  mutate(inner_rec_phi = round(inner_rec_phi,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))



CMR_juv %>% 
  spread_draws(gamma, gamma_rec) %>% 
  median_qi(outer_rec_gamma = gamma + gamma_rec, .width = 0.9) %>% 
  mutate(outer_rec_gamma = round(outer_rec_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))


CMR_juv %>% 
  spread_draws(gamma, gamma_rec,gamma_rec_interac_inout) %>% 
  median_qi(inner_rec_gamma = gamma + gamma_rec +  gamma_rec_interac_inout, .width = 0.9) %>% 
  mutate(inner_rec_gamma = round(inner_rec_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(gamma, gamma_rec,gamma_rec_interac_inout) %>% 
  mutate(outer_rec_gamma = gamma + gamma_rec) %>% 
  mutate(inner_rec_gamma = gamma + gamma_rec +  gamma_rec_interac_inout) %>% 
  mutate(prop_rec_g = outer_rec_gamma/inner_rec_gamma) %>% 
  median_qi(prop_rec_g, .width =0.9)
launch_shinystan(CMR_juv)

CMR_juv %>% 
  spread_draws(gamma_rec_interac_inout) %>% 
  median_qi(inner_rec_gamma_diff = gamma_rec_interac_inout, .width = 0.9) %>% 
  mutate(inner_rec_gamma_diff = round(inner_rec_gamma_diff,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_juv %>% 
  spread_draws(gamma, gamma_rec,gamma_rec_interac_inout) %>% 
  median_qi(metapop_rec_gamma = gamma + gamma_rec +  gamma_rec_interac_inout*0.5, .width = 0.9) %>% 
  mutate(metapop_rec_gamma = round(metapop_rec_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

#survival prob from nestling to recruit

CMR_juv %>% 
  spread_draws(mu_phi,phi_fledge_interac_inout, mu_rec_phi,phi_rec_interac_inout) %>% 
  mutate(outer_nest_to_rec_phi = invlogit(mu_phi) * invlogit(mu_phi +mu_rec_phi)) %>% 
  mutate(inner_fledge_phi = invlogit(mu_phi + phi_fledge_interac_inout)) %>% 
  mutate(inner_rec_phi = invlogit(mu_phi+mu_rec_phi + phi_rec_interac_inout)) %>% 
  mutate(inner_nest_to_rec_phi = inner_fledge_phi * inner_rec_phi) %>% 
  mutate(diff_outin =  outer_nest_to_rec_phi - inner_nest_to_rec_phi) %>% 
  median_qi(outer_nest_to_rec_phi,inner_nest_to_rec_phi,diff_outin, .width = 0.9)



#adults
CMR_adult_survival %>% 
  get_variables() %>% 
  head(n=12)



CMR_adult_survival %>% 
  spread_draws(mu_phi) %>% 
  median_qi(outer_phi = invlogit(mu_phi) , .width = 0.9) %>% 
  mutate(outer_phi = round(outer_phi,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_adult_survival %>% 
  spread_draws(mu_phi, B_in_out) %>% 
  median_qi(inner_phi = invlogit(mu_phi + B_in_out) , .width = 0.9) %>% 
  mutate(inner_phi = round(inner_phi,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

#Calculate survival probability of the average sex at the average age. I've double checked that
#multiplying with 0.5 works. It lowers the male part of the intercept by half.
mean_age<-mean(data_observed$Least_age-1) #-1 so that recruits are 0. That is how we modelled age.

CMR_adult_survival %>% 
  get_variables() %>% 
  head(n=12)

(avg_outer<-CMR_adult_survival %>% 
    spread_draws(mu_phi,B_sex,B_age, B_age_sq,interac_age_sex) %>% 
    median_qi(outer_average_phi = 
                
                invlogit(mu_phi + B_sex*0.5 + B_age * mean_age + B_age_sq * mean_age^2+ interac_age_sex * mean_age *0.5) 
              
              , .width = 0.9) %>% 
    mutate(outer_average_phi = round(outer_average_phi,2),
           .lower= round(.lower,2), 
           .upper = round(.upper,2)))

(avg_inner<-CMR_adult_survival %>% 
    spread_draws(mu_phi,B_sex,B_age, B_age_sq,interac_age_sex, B_in_out) %>% 
    median_qi(inner_average_phi = 
                
                invlogit(mu_phi + B_sex*0.5 + B_age * mean_age + B_age_sq * mean_age^2+ interac_age_sex * mean_age *0.5 + B_in_out) 
              
              , .width = 0.9) %>% 
    mutate(inner_average_phi = round(inner_average_phi,2),
           .lower= round(.lower,2), 
           .upper = round(.upper,2)))

CMR_adult_survival %>% 
  spread_draws(mu_phi,B_sex,B_age, B_age_sq,interac_age_sex, B_in_out) %>% 
  mutate(inner_average_phi = 
           invlogit(mu_phi + B_sex*0.5 + B_age * mean_age + B_age_sq * mean_age^2+ interac_age_sex * mean_age *0.5 + B_in_out)) %>% 
  mutate(outer_average_phi = 
           invlogit(mu_phi + B_sex*0.5 + B_age * mean_age + B_age_sq * mean_age^2+ interac_age_sex * mean_age *0.5)) %>% 
  median_qi(prop_surv  =inner_average_phi/outer_average_phi, .width = 0.9) %>% 
  mutate(prop_surv = round(prop_surv,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_adult_survival %>% 
  spread_draws(gamma) %>% 
  median_qi(outer_gamma = gamma , .width = 0.9) %>% 
  mutate(outer_gamma = round(outer_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_adult_survival %>% 
  spread_draws(gamma, interac_in_out_gamma) %>% 
  median_qi(inner_gamma = gamma +interac_in_out_gamma , .width = 0.9) %>% 
  mutate(inner_gamma = round(inner_gamma,2),
         .lower= round(.lower,2), 
         .upper = round(.upper,2))

CMR_adult_survival %>% 
  spread_draws(gamma, interac_in_out_gamma) %>%
  mutate(g_meta = gamma + interac_in_out_gamma*0.5) %>% 
  median_qi(round(g_meta,2), .width = 0.9)


#Descriptive statistics ####

##More predictions adult fitness ####
bernt_zNR_nb_brm_noindslope %>% 
  get_variables() %>% 
  head(n=12)
bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept, b_Least_age, b_Least_age_sq, b_in_out, b_sex, `b_Least_age:sex`) %>% 
  mutate(avg_fit =b_Intercept + b_Least_age*mean_age + b_Least_age_sq*mean_age^2 + b_in_out) %>% 
  median_qi(exp(avg_fit-log(2)), .width = 0.9)


bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept) %>% 
  median_qi(exp(b_Intercept-log(2)), .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>%
  mutate(inner_g = b_zNR + `b_zNR:in_out`) %>% 
  mutate(prop_g = b_zNR/inner_g) %>% 
  median_qi(prop_g, .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR) %>%
  mutate(outer_g = b_zNR) %>% 
  median_qi(outer_g, .width = 0.9)



bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>%
  mutate(inner_g = b_zNR + `b_zNR:in_out`) %>% 
  median_qi(inner_g, .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept, b_zNR) %>% 
  mutate(pred =b_Intercept + b_zNR) %>% 
  median_qi(exp(pred-log(2)), .width = 0.9)



bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR, `b_zNR:in_out`) %>%
  mutate(meta_g = b_zNR + `b_zNR:in_out`*0.5) %>% 
  mutate(fitness_1zNR = b_Intercept + b_in_out*0.5 + meta_g*1) %>% 
  median_qi(exp(fitness_1zNR-log(2)), .width = 0.9)



bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR, `b_zNR:in_out`) %>%
  mutate(meta_g = b_zNR + `b_zNR:in_out`*0.5) %>% 
  mutate(fitness_equil = b_Intercept + b_in_out*0.5) %>% 
  mutate(fitness_1zNR = fitness_equil + meta_g*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(prop, .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR) %>%
  mutate(fitness_equil = b_Intercept) %>% 
  mutate(fitness_1zNR = fitness_equil + b_zNR*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(prop, .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_in_out, b_zNR, `b_zNR:in_out`) %>%
  mutate(inner_g = b_zNR + `b_zNR:in_out`) %>% 
  mutate(fitness_equil = b_Intercept + b_in_out) %>% 
  mutate(fitness_1zNR = fitness_equil + inner_g*1) %>% 
  mutate(prop = exp(fitness_1zNR-log(2))/exp(fitness_equil-log(2))) %>% 
  median_qi(prop, .width = 0.9)


bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_in_out) %>%
  mutate(mean_fitness = b_Intercept + b_in_out*0.5) %>% 
  median_qi(exp(mean_fitness-log(2)), .width = 0.9)

rec_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>%
  mutate(g_metapop = b_zNR + `b_zNR:in_out`*0.5) %>% 
  median_qi(g_metapop, .width = 0.9)

rec_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>%
  mutate(g_inner = b_zNR + `b_zNR:in_out`) %>% 
  median_qi(g_inner, .width = 0.9)


rec_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR) %>%
  mutate(g_outer = b_zNR) %>% 
  median_qi(g_outer, .width = 0.9)

#percent difference in gamma of recruit prod between inner and outer width credible intervals 
rec_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>%
  mutate(g_inner = b_zNR + `b_zNR:in_out`) %>% 
  mutate(prop_diff = b_zNR/g_inner) %>% 
  median_qi(prop_diff, .width = 0.9) %>% 
  mutate(prop_diff  = (prop_diff -1)*100) %>% 
  mutate(.lower = (.lower-1)*100) %>% 
  mutate(.upper = (.upper-1)*100)


rec_zNR_nb_brm_noindslope %>% 
  spread_draws(b_Intercept,b_Least_age, `b_Least_age_sq`, b_in_out) %>%
  mutate(rec_inner = b_Intercept + b_Least_age*mean_age +b_Least_age_sq*mean_age^2) %>% 
  median_qi(exp(rec_inner-log(2)), .width = 0.9)



##recapture probs for each island and year ####
recapFY<-CMR_adult_survival %>% 
  spread_draws(p[Location,Year]) %>% 
  mutate(Year = Year+ 1998) %>% 
  group_by(Location,Year) %>% 
  median_qi(p,.width=0.9) %>% 
  ungroup() %>% 
  mutate(p = round(p,2)) %>% 
  mutate(.lower = paste0("(", round(.lower,2), ",")) %>% 
  mutate(.upper = paste0(round(.upper,2), ")")) %>% 
  unite(p, p,.lower, .upper, sep = " ") %>% 
  dplyr::select(Location, Year, p) %>%
  pivot_wider(names_from = Year, values_from = p)

isl_names<-data_observed %>% 
  dplyr::select(Location) %>% 
  arrange(Location) %>% 
  unique()

recapFY$Location <- isl_names$Location

View(recapFY)

recapF<-cjs_temp_raneff2 %>% 
  spread_draws(p[Location, Year]) %>% 
  group_by(Location) %>% 
  median_qi(p,.width=0.9) %>% 
  ungroup() %>% 
  mutate(p = round(p,2)) %>% 
  mutate(.lower = paste0("(", round(.lower,2), ",")) %>% 
  mutate(.upper = paste0(round(.upper,2), ")")) %>% 
  unite(p, p,.lower, .upper, sep = " ") %>% 
  dplyr::select(Location, p)

recapF$Location <- isl_names$Location
recapF

recapY<-cjs_temp_raneff2 %>% 
  spread_draws(p[Location, Year]) %>% 
  mutate(Year = Year + 1998) %>% 
  group_by(Year) %>% 
  median_qi(p,.width=0.9) %>% 
  ungroup() %>% 
  mutate(p = round(p,2)) %>% 
  mutate(.lower = paste0("(", round(.lower,2), ",")) %>% 
  mutate(.upper = paste0(round(.upper,2), ")")) %>% 
  unite(p, p,.lower, .upper, sep = " ") %>% 
  dplyr::select(Year, p)


recapY

## Variance in gamma explained by inner outer ####

bernt_zNR_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_zNR_nb_brm_noindslope.rds"))
bernt_zNR_nb_brm_no_in_out_slope<-readRDS("Workspace backup/bernt_zNR_nb_brm_no_in_out_slope.rds")
bernt_n_nb_brm_no_in_out_slope<-readRDS("Workspace backup/bernt_n_nb_brm_no_in_out_slope.rds")
bernt_n_nb_brm_noindslope<-readRDS(file("Workspace backup/bernt_n_nb_brm_noindslope.rds"))

zNR_no_inout<-brm2table_short(bernt_zNR_nb_brm_no_in_out_slope, type = "median", width = 0.9)
zNR_with_inout<-brm2table_short(bernt_zNR_nb_brm_noindslope, type = "median", width = 0.9)
n_no_inout<-brm2table_short(bernt_n_nb_brm_no_in_out_slope, type = "median", width = 0.9)
n_with_inout<-brm2table_short(bernt_n_nb_brm_noindslope, type ="median" , width =0.9)

zNR_no_inout %>% 
  filter(Parameter == "Sigma2_Isl_gamma")
zNR_with_inout%>% 
  filter(Parameter == "Sigma2_Isl_gamma")
n_no_inout%>% 
  filter(Parameter == "Sigma2_Isl_gamma")
n_with_inout%>% 
  filter(Parameter == "Sigma2_Isl_gamma")


bernt_zNR_nb_brm_no_in_out_slope %>% 
  spread_draws(sd_Location__zNR) %>% 
  median_qi(sd_Location__zNR^2, .width = 0.9)

bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(sd_Location__zNR) %>% 
  median_qi(sd_Location__zNR^2, .width = 0.9)

sd_no_inout<-bernt_zNR_nb_brm_no_in_out_slope %>% 
  spread_draws(sd_Location__zNR)

sd_with_inout<-bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(sd_Location__zNR)
par(mfrow = c(1,2))
hist((sd_with_inout$sd_Location__zNR)^2)
hist((sd_no_inout$sd_Location__zNR)^2)

prop<-(sd_with_inout$sd_Location__zNR^2)/(sd_no_inout$sd_Location__zNR^2)

median_qi(prop, .width= 0.9)
##Total variance and proportion of variance explained by gamma ####
#total variance in log mean fitness
total_var_mean_r2s<-data_observed %>% 
  group_by(Year_Location) %>% 
  summarise(log_mean_r2s = log(mean(r2s))) %>% 
  summarise(var_log_mean_r2s = var(log_mean_r2s))
mean(data_observed$zNR)

#total variance in log mean fitness due to gamma
var_zNR_outer<-data_observed %>% 
  filter(in_out == 0) %>% 
  dplyr::select(Year_Location, zNR) %>% 
  unique() %>% 
  summarise(var_zNR_outer = mean(zNR))

var_zNR_inner<-data_observed %>% 
  filter(in_out == 1) %>% 
  dplyr::select(Year_Location, zNR) %>% 
  unique() %>% 
  summarise(var_zNR_inner = var(zNR))

var_zNR_tot<-data_observed %>% 
  dplyr::select(Year_Location, zNR) %>% 
  unique() %>% 
  summarise(var_zNR_tot = var(zNR))

var_slopes = as.data.frame(VarCorr(bernt_zNR_nb_brm_noindslope,summary = FALSE)$Location$cov)[,4]

total_gamma_var<-bernt_zNR_nb_brm_noindslope %>% 
  spread_draws(b_zNR, `b_zNR:in_out`) %>% 
  mutate(outer_gamma_var = b_zNR^2 *var_zNR_outer$var_zNR_outer) %>% 
  mutate(inner_gamma_var = (b_zNR +`b_zNR:in_out`)^2 * var_zNR_inner$var_zNR_inner) %>% 
  mutate(random_slopes_var = var_slopes * var_zNR_tot$var_zNR_tot) %>%
  #median_qi(prop_inout_gamma =(outer_gamma_var + inner_gamma_var)/ (outer_gamma_var + inner_gamma_var + random_slopes_var), .width = 0.9) %>% #can't be right
  mutate(total_gamma_var = outer_gamma_var + inner_gamma_var + random_slopes_var) %>% 
  median_qi(total_gamma_var/total_var_mean_r2s$var_log_mean_r2s, .width = 0.9)

##Single islands correlation between stochasticty and r0 and g ####

bernt_zmc_N_single_island_mods<-readRDS(file("Workspace backup/bernt_zmc_N_single_island_mods.rds"))

island_vars<-as.data.frame(matrix(NA, 11, 6))
for(i in 1:11){
  island_vars[i,]<-median_qi(VarCorr(bernt_zmc_N_single_island_mods[[i]], summary = FALSE)$Year$sd^2)
  island_vars[i,7] <- names(bernt_zmc_N_single_island_mods)[i]
}
colnames(island_vars) <- c("Sigma2_Year", ".lower", ".upper", ".width", ".point", ".interval","Location")
island_vars[,1:3] <- round(island_vars[,1:3],3)
island_vars %>% 
  arrange(Sigma2_Year)
island_vars



divide_intercepts = "yes"

post_frame= get_single_island_estimates(bernt_zmc_N_single_island_mods,divide_intercepts, summary = TRUE)
post_frame

g_var<-island_vars %>% 
  left_join(post_frame, join_by("Location" == "Island")) %>% 
  dplyr::select(Location, Sigma2_Year, g)

plot(g ~ Sigma2_Year, data = g_var)

g_var %>% 
  arrange(Sigma2_Year)

island_N_means<-data_observed %>% 
  dplyr::select(Year, Location, estimated_pop_size) %>% 
  unique() %>% 
  group_by(Location) %>% 
  summarise(Mean_N_isl = mean(estimated_pop_size)) %>% 
  ungroup()
sd_N<-data_observed %>% 
  dplyr::select(Location, Year, N) %>% 
  unique() %>% 
  summarise(sd(N))

zmc_single_backtransformed_r0 <- post_frame %>% 
  left_join(island_N_means, by = c("Island" = "Location")) %>% 
  mutate(r0_b = r0-(g*Mean_N_isl/ sd_N$`sd(N)`)) %>% 
  dplyr::select(Island,r0_b)

g_var

r0_g_var <- g_var %>% 
  left_join(zmc_single_backtransformed_r0, by = c("Location" = "Island"))

r0_g_var %>% 
  arrange(Sigma2_Year)

plot(r0_b ~ Sigma2_Year, data = r0_g_var) 



##Variance in r0 and g in models n, N_m, single N_m ####

density_variable = "n"
fitness = "r2s"
divide_intercepts= "yes"
average_metapars ="yes"

post_nb_full_n<-brm_blups_islands(bernt_n_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)
post_nb_full_n

density_variable = "zmc_N"
fitness = "r2s"
divide_intercepts = "yes"
average_metapars ="yes"
mean_centred = "yes"


post_nb_full_Nm<-brm_blups_islands(bernt_zNR_nb_brm_noindslope, density_variable, fitness, divide_intercepts,average_metapars)
post_nb_full_Nm

bernt_zNR_single_island_mods<-readRDS(file("Workspace backup/bernt_zNR_single_island_mods.rds"))
bernt_zNR_single_island_mods[[1]]
divide_intercepts = "yes"

post_frame_single= get_single_island_estimates(bernt_zNR_single_island_mods,divide_intercepts, summary = TRUE)
post_frame_single

par(mfrow = c(3,1))
par(mar = c(4,4,4,4))




post_nometa_n<-post_nb_full_n %>% 
  filter(Island !="meta-population")
var(exp(post_nometa_n$r0))
#This is where we calculate the variance in slopes presented in the text (interaction of inner_outer and g was added to g in the brm_blups_islands function)
var(exp(post_nometa_n$g))
exp(post_nometa_n$r0)
max(exp(post_nometa_n$r0))



post_nometa_n <- post_nometa_n %>% 
  left_join(island_N_means, by =c("Island" = "Location"))

plot(r0  ~ Mean_N_isl, data = post_nometa_n)

island_N_means<-data_observed %>% 
  dplyr::select(Location, Mean_N_isl) %>% 
  unique()

r0_backtransformed_zmc <- post_nb_full_Nm %>% 
  left_join(island_N_means, by = c("Island"= "Location")) %>% 
  mutate(r0 = r0 - g*Mean_N_isl/40) %>% #the "old" r0 is actually r_mu
  mutate(r0 = r0.lower - g*Mean_N_isl/40) %>% 
  mutate(r0 = r0.upper - g*Mean_N_isl/40) %>% 
  filter(Island != "meta-population")
r0_backtransformed_zmc %>% 
  mutate(exp(r0))

var(r0_backtransformed_zmc$r0)
max(exp(r0_backtransformed_zmc$r0))

plot(r0  ~ Mean_N_isl, data = r0_backtransformed_zmc)

#This is where we calculate the variance in slopes presented in the text (interaction of inner_outer and g was added to g in the brm_blups_islands function)
post_nb_full_Nm %>% 
  filter(Island != "meta-population") %>% 
  dplyr::select(g) %>% 
  var()
#single
zmc_single_backtransformed <- post_frame_single %>% 
  left_join(island_N_means, by = c("Island" = "Location")) %>% 
  mutate(r0 = r0-(g*Mean_N_isl/ 40)) %>% 
  mutate(r0.lower = r0.lower - g*Mean_N_isl/40) %>% 
  mutate(r0.lower = r0.upper - g*Mean_N_isl/40)

zmc_single_backtransformed %>% 
  mutate(exp(r0))


var(zmc_single_backtransformed$r0)
var(r0_backtransformed_zmc$r0)
var(post_nometa_n$r0)
var(post_frame$g)

zmc_single_backtransformed %>% 
  mutate(exp(r0))
max(exp(zmc_single_backtransformed$r0))
var(exp(zmc_single_backtransformed$r0))

plot(r0 ~ Mean_N_isl, data = zmc_single_backtransformed)

post_nb_full_Nm %>% 
  filter(Island != "meta-population") %>% 
  dplyr::select(g) %>% 
  var()
var(zmc_single_backtransformed$g)
var(post_nometa_n$g) #exponentiate or not?

##r0 tables islands ####


tables <- list()

tables[[1]] <- r0_backtransformed_zmc
tables[[2]] <- zmc_single_backtransformed
tables[[3]] <- post_nometa_n
varextr <- function(x){var(x[,5])}
varextr(r0_backtransformed_zmc)
lapply(tables, varextr)

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

##Fitness at the mean density tables islands ####

post_nb_full_Nm_nometa<-post_nb_full_Nm %>% 
  filter(Island != "meta-population")

post_nometa_n_meanpopsize<-post_nometa_n %>% 
  mutate(r0 = r0 + g*log(Mean_N_isl)) %>% 
  mutate(r0.lower = r0.lower + g*log(Mean_N_isl)) %>% 
  mutate(r0.upper = r0.upper + g*log(Mean_N_isl))

tables <- list()

tables[[1]] <- post_nb_full_Nm_nometa
tables[[2]] <- post_frame_single
tables[[3]] <- post_nometa_n_meanpopsize
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



##  Correlation between pop size and gamma and r0 ####
r0_backtransformed_zmc
zmc_single_backtransformed
post_nometa_n

inout<-c("aldra", "gjerøy", "hestmannøy", "indre kvarøy", "lurøy-onøy", "nesøy")

gr0_zNr_meta<-r0_backtransformed_zmc %>% 
  mutate(in_out = ifelse(Island %in% inout, "inner", "outer"))


gr0_zNr_single<-zmc_single_backtransformed %>% 
  mutate(in_out = ifelse(Island %in% inout, "inner", "outer"))

gr0n_meta<-post_nometa_n %>% 
  mutate(in_out = ifelse(Island %in% inout, "inner", "outer"))



plot.new()
par(mfrow = c(3,3))
par(mar= c(2,2,2,2))
par(oma= c(6,1,4,2))

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression('Logistic, meta-population model (N'[m]*')'), cex = 1.4)

plot(g~Mean_N_isl, data = gr0_zNr_meta, xlab="", ylab="",cex = 1.5,
     col=as.factor(gr0_zNr_meta$in_out),pch= as.numeric(as.factor(gr0_zNr_meta$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0_zNr_meta$in_out)), col = 1:2, cex = 0.1, pch=unique(as.numeric(as.factor(gr0_zNr_meta$in_out))))
plot(r0~Mean_N_isl, data = gr0_zNr_meta, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_meta$in_out),pch= as.numeric(as.factor(gr0_zNr_meta$in_out)), cex.axis=1.3)

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression('Logistic, single population model (N'[m]*')'), cex = 1.4)

plot(g~Mean_N_isl, data = gr0_zNr_single, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_single$in_out), pch= as.numeric(as.factor(gr0_zNr_single$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0_zNr_single$in_out)), col = 1:2, cex = 1, pch=unique(as.numeric(as.factor(gr0_zNr_single$in_out))))
plot(r0~Mean_N_isl, data = gr0_zNr_single, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0_zNr_single$in_out), pch= as.numeric(as.factor(gr0_zNr_single$in_out)), cex.axis=1.3)

plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
     xaxt = "n", yaxt = "n")
text(x = 5,y = 5,expression("Gompertz meta-population model (n)"), cex = 1.4)


plot(g~Mean_N_isl, data = gr0n_meta, xlab="", ylab="", cex = 1.5,
     col=as.factor(gr0n_meta$in_out),pch= as.numeric(as.factor(gr0n_meta$in_out)), cex.axis=1.3)
#legend('bottomright', legend = levels(as.factor(gr0n_meta$in_out)),  col = 1:2, cex = 0.1, pch=unique(as.numeric(as.factor(gr0n_meta$in_out))))
plot(r0~Mean_N_isl, data = gr0n_meta, xlab="", ylab="",cex = 1.5,
     col=as.factor(gr0n_meta$in_out),pch= as.numeric(as.factor(gr0n_meta$in_out)), cex.axis=1.3)
legend("bottomright",
       legend = levels(as.factor(gr0n_meta$in_out)), 
       col = 1:2, 
       xpd = TRUE, cex = 1.2,pch=unique(as.numeric(as.factor(gr0n_meta$in_out))), seg.len=1)

mtext(expression(gamma), side = 3, outer = T, cex = 1.3)
mtext(expression("r"[0]), side = 3, outer = T, adj = 0.85, cex = 1.3)
mtext("Island-specific mean density", side = 1, outer = T, adj = 0.7, cex = 1.3,line = 2)


##Inner outer mean pop size ####
data_observed$N <- data_observed$estimated_pop_size
popsizes<-data_observed %>% 
  dplyr::select(Location, Year, N, in_out) %>% 
  unique()
popsizes %>% 
  group_by(in_out) %>% 
  summarise(mean(N))
popsizes %>% 
  group_by(in_out,Location) %>% 
  summarise(mean(N))
popsizes %>% 
  group_by(Location,in_out) %>% 
  summarise(var(N))



hist(popsizes$N)
in_out_means<-glmer.nb(N~in_out + (1|Location), data= popsizes)
summary(in_out_means)

conf_inout<-confint(in_out_means)
#inner pop size
exp(fixef(in_out_means)[1] + VarCorr(in_out_means)$Location[1]/2)
exp(conf_inout[2,1] + VarCorr(in_out_means)$Location[1]/2)
exp(conf_inout[2,2] + VarCorr(in_out_means)$Location[1]/2)
#outer pop size

exp(fixef(in_out_means)[1] + fixef(in_out_means)[2]+  + VarCorr(in_out_means)$Location[1]/2)
exp(conf_inout[2,1] + conf_inout[3,1]+ VarCorr(in_out_means)$Location[1]/2)
exp(conf_inout[2,2] + conf_inout[3,2] +VarCorr(in_out_means)$Location[1]/2)

var(subset)
##Inner outer mean gamma and mean r0 single island models ####
bernt_zNR_single_island_mods<-readRDS(file("Workspace backup/bernt_zNR_single_island_mods.rds"))
divide_intercepts="yes"

post_frame_post= get_single_island_estimates(bernt_zNR_single_island_mods,divide_intercepts, summary = FALSE)

single_inout<-post_frame_post %>% 
  mutate(in_out = ifelse(Island %in% c("gjerøy", "hestmannøy", "indre kvarøy", "aldra", "lurøy-onøy", "nesøy"), "inner", "outer"))

single_inout %>%   
  group_by(in_out) %>% 
  median_qi(r0,g)

head(single_inout)
single_inout %>% 
  ggplot(aes(y = Island, x = r0)) +
  stat_halfeye(.width = c(0.95, 1), point_interval="median_qi",aes(fill = after_stat(level)))
single_inout %>% 
  ggplot(aes(y = in_out, x = r0)) +
  stat_halfeye(.width = c(0.95, 1), point_interval="median_qi",aes(fill = after_stat(level)))+
  theme_bw()

## Population size per island over time figure ####
par(mfrow = c(1,1))
isl_pop<-data_observed %>% 
  dplyr::select(Location, Year, estimated_pop_size) %>% 
  unique() %>% 
  arrange(Location, Year)
plot(1,type ="n", ann = FALSE,
     ylim = c(0,max(data_observed$estimated_pop_size)),
     xlim = c(min(isl_pop$Year), max(isl_pop$Year)), xaxt = "n", ylab = "N")

yrs<-seq(1998, 2019, by = 3)

axis(1, at=yrs)
title(xlab = "Year",ylab = "N")
island_colours<-c(rainbow(8),"#00B3FF","#FF004D","#0066FF")

for(i in 1:unique(length(isl_pop$Location))){
  single_island <- subset(isl_pop, Location == unique(isl_pop$Location)[i])
  points(single_island$estimated_pop_size ~ single_island$Year, type = "l",lty = i, col = island_colours[i], lwd = 2)
}

legend("topright", legend = unique(isl_pop$Location),
       
       col=c(island_colours),lty =c(1:11), lwd = 2,
       title="Island",  bg='grey')

## Correlation between estimated pop size (pre-breeding census) and our fitness estimates (r2s/2) ####
fits<-data_observed %>% 
  group_by(Location, Year) %>% 
  mutate(tot_fitness = sum(r2s)/2) %>% 
  ungroup %>% 
  dplyr::select(Location, Year,tot_fitness) %>% 
  mutate(Year = Year +1 ) %>% 
  filter(Year != 2020) %>% 
  unique()

pops <- data_observed %>% 
  dplyr::select(Location, Year, estimated_pop_size)

pops <- data_observed %>% 
  dplyr::select(Location, Year, estimated_pop_size) %>% 
  #filter(Year != 1998) %>% 
  unique()
fits<-fits %>% 
  left_join(pops, by = c("Location", "Year")) %>% 
  arrange(Location, Year)
par(mfrow=c(6,2))
par(mar = c(1,2,3,1))
par(oma = c(4,4,1,1))
for(i in 1:11){
  single_island <- subset(fits, Location == unique(fits$Location)[i])
  
  plot(tot_fitness ~ estimated_pop_size, data = single_island, main = unique(single_island$Location), xlab = "", ylab ="")
}

yrs<-seq(1999, 2019, by = 4)

plot.new()
for(i in 1:11){
  single_island <- subset(fits, Location == unique(fits$Location)[i])
  
  plot(estimated_pop_size~Year, data = single_island, main = unique(single_island$Location), xlab = "", ylab ="", type ="l", xaxt = "n")
  axis(1, at=yrs)
  points(tot_fitness~Year, data = single_island, col = "red", type = "l")
} 

plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', xlab="", ylab="")
legend('center',
       legend = c("Pre-breeding census", "Recruitment + survival"), 
       col = c("black","red"), 
       lwd=2, xpd = TRUE, cex = 2, seg.len=1, bty = 'n')
mtext(text="Year",side=1,line=2,outer=TRUE)
mtext(text="Number of adults",side=2,line=1,outer=TRUE)

