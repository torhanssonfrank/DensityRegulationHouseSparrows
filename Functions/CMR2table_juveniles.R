#CMR to table juveniles

#I want one table for nestling to fledgling and one table for fledgling to recruit
# so I need to add the intercepts together to present the actual stage-specific
#estimates rather than the differences. But the downside is that you cannot see
#if the differences are significant. But maybe I just put the raw output in the supmat.

CMR2table_juveniles<-function(model,paramlatex,reverse_inner_outer,add_difference_to_intercept, width){

  library(tidyverse)
  library(tidybayes)
  
CMR_juv_post<-as.data.frame(model)

CMR_juv_post<-CMR_juv_post %>% 
  dplyr::select(any_of(paramlatex$Parameter))

CMR_juv_post<- CMR_juv_post %>% 
  rename(mu_p_fledge = mu_p) %>% 
  rename(mu_fledge_phi = mu_phi) %>% 
  rename(gamma_fledge = gamma)


if(reverse_inner_outer == "yes"){
  CMR_juv_post<-CMR_juv_post %>% 
    mutate(mu_fledge_phi_new = mu_fledge_phi + phi_fledge_interac_inout) %>% 
    mutate(phi_fledge_interac_inout_new = mu_fledge_phi - mu_fledge_phi_new) %>% 
    mutate(gamma_fledge_new = gamma_fledge + gamma_fledge_interac_inout) %>% 
    mutate(gamma_fledge_interac_inout_new = gamma_fledge - gamma_fledge_new) %>% 
    
    mutate(mu_rec_phi_outer_full = mu_fledge_phi + mu_rec_phi + phi_rec_interac_inout) %>% #maybe something wrong here
    mutate(mu_rec_phi_inner_full = mu_fledge_phi + mu_rec_phi) %>% 
    
    mutate(mu_rec_phi_new = mu_rec_phi_outer_full - mu_fledge_phi_new) %>% #should be difference between outer recs and outer fledglings
    mutate(phi_rec_interac_inout_new = mu_rec_phi_inner_full - mu_rec_phi_outer_full) %>% #should be difference to outer recruits
  
    mutate(gamma_rec_outer_full = gamma_fledge + gamma_rec + gamma_rec_interac_inout) %>% 
    mutate(gamma_rec_inner_full = gamma_fledge + gamma_rec) %>% 
    
    
    mutate(gamma_rec_new = gamma_rec_outer_full - gamma_fledge_new) %>% #should be difference between outer recs and outer fledglings
    mutate(gamma_rec_interac_inout_new = gamma_rec_inner_full - gamma_rec_outer_full) %>% #should be difference to outer recruits
    
    dplyr::select(-mu_fledge_phi,-phi_fledge_interac_inout,-gamma_fledge,-gamma_fledge_interac_inout,
                  -mu_rec_phi, -phi_rec_interac_inout,-gamma_rec,-gamma_rec_interac_inout, 
                  -mu_rec_phi_outer_full, - mu_rec_phi_inner_full, -gamma_rec_outer_full, -gamma_rec_inner_full) %>% 
    rename(
      mu_fledge_phi= mu_fledge_phi_new,
      phi_fledge_interac_inout = phi_fledge_interac_inout_new,
      gamma_fledge= gamma_fledge_new,
      gamma_fledge_interac_inout=gamma_fledge_interac_inout_new,
      
      mu_rec_phi = mu_rec_phi_new,
      phi_rec_interac_inout = phi_rec_interac_inout_new,
      gamma_rec = gamma_rec_new,
      gamma_rec_interac_inout = gamma_rec_interac_inout_new)
  
} else if(reverse_inner_outer == "no"){
  CMR_juv_post <- CMR_juv_post
  
} else{
  print("Error, reverse_inner_outer must be yes or no!")
}


#calculate the actual point estimates for recruits (add fledge estimate to difference)

if(add_difference_to_intercept == "yes"){
CMR_juv_res<-CMR_juv_post %>% 
  mutate(mu_rec_phi = mu_rec_phi + mu_fledge_phi) %>% 
  mutate(mu_rec_p = mu_rec_p + mu_p_fledge)

} else if(add_difference_to_intercept == "no"){
  CMR_juv_res<-CMR_juv_post
  
} else{
  print("Error, add_difference_to_intercept must be yes or no!")
}

  
  
CMR_juv_point_estimates<-CMR_juv_res %>% 
  pivot_longer(1:ncol(CMR_juv_res), names_to = "Parameter", values_to = "Estimate") %>% 
  group_by(Parameter) %>% 
  median_qi(Estimate, .width = width)



fledglings <- CMR_juv_point_estimates %>% 
  filter(str_detect(Parameter, "fledge"))

recruits <- CMR_juv_point_estimates %>% 
  filter(str_detect(Parameter, "rec"))

juv_tables <- list()

juv_tables[[1]] <- fledglings
juv_tables[[2]] <- recruits

juv_tables

}

