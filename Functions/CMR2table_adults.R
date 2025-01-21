#Create table from CMR adult survival density regulation
library(tidyverse)
library(tidybayes)

CMR2table_adults<-function(model,paramlatex,reverse_inner_outer, width){
CMR_adult_posterior<-as.data.frame(model)

CMR_adult_posterior<-CMR_adult_posterior %>% 
  dplyr::select(any_of(paramlatex$Parameter))


if(reverse_inner_outer == "yes"){
  CMR_adult_posterior<-CMR_adult_posterior %>% 
    mutate(mu_phi_new = mu_phi + B_in_out) %>% 
    mutate(B_in_out_new = mu_phi - mu_phi_new) %>% 
    mutate(gamma_new = gamma + interac_in_out_gamma) %>% 
    mutate(interac_in_out_gamma_new = gamma - gamma_new) %>% 
    dplyr::select(-mu_phi,-B_in_out,-interac_in_out_gamma,-gamma) %>% 
    rename(
      mu_phi= mu_phi_new,
      B_in_out = B_in_out_new,
      gamma= gamma_new,
      interac_in_out_gamma=interac_in_out_gamma_new)
} else if(reverse_inner_outer == "no"){
  CMR_adult_posterior <- CMR_adult_posterior
} else{
  print("Error, reverse_inner_outer must be yes or no!")
}




CMR_adult_posterior_long<-CMR_adult_posterior %>% 
  pivot_longer(1:ncol(CMR_adult_posterior), names_to = "Parameter", values_to = "Estimate")



CMR_adult_table<-CMR_adult_posterior_long %>% 
  group_by(Parameter) %>% 
  median_qi(Estimate, .width = width)
CMR_adult_table
}