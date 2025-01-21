#Functions to generate tables and figures from stan and brms outputs

table_wrangler<-function(reproduction, survival, w2s, paramlatex){
  library(Matrix)
  library(tidyverse)
  
  colnames(reproduction)[colnames(reproduction) == ".lower"] <-paste0(print((1-reproduction$.width[1])/2*100),"%")
  colnames(reproduction)[colnames(reproduction) == ".upper"] <-paste0(print(100-(1-(reproduction$.width[1]))/2*100),"%")
  
  
  reproduction<-reproduction %>% 
    mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>% 
    mutate(Reproduction = paste0(Estimate, " (", `2.5%`,", ", `97.5%`,")")) %>% 
    dplyr::select(Parameter, Reproduction)
  
  colnames(survival)[colnames(survival) == ".lower"] <-paste0(print((1-survival$.width[1])/2*100),"%")
  colnames(survival)[colnames(survival) == ".upper"] <-paste0(print(100-(1-(survival$.width[1]))/2*100),"%")
  
  
  survival<-survival %>% 
    mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>%
    mutate(Survival = paste0(Estimate, " (", `2.5%`,", ", `97.5%`,")")) %>% 
    dplyr::select(Parameter, Survival)
  
  colnames(w2s)[colnames(w2s) == ".lower"] <-paste0(print((1-w2s$.width[1])/2*100),"%")
  colnames(w2s)[colnames(w2s) == ".upper"] <-paste0(print(100-(1-(w2s$.width[1]))/2*100),"%")
  
  
  w2s<-w2s%>% 
    mutate(across(where(is.numeric), ~ round(.x, digits = 2))) %>% 
    mutate(`Reproduction + Survival` = paste0(Estimate, " (", `2.5%`,", ", `97.5%`,")")) %>% 
    dplyr::select(Parameter, `Reproduction + Survival`)
  
  fit_table<-reproduction %>% 
    left_join(survival) %>% 
    left_join(w2s) 
  
  
  
  
  fit_table<-fit_table %>% 
    left_join(paramlatex) %>% 
    dplyr::select(-Parameter) %>% 
    rename(Parameter = Parameter.latex)%>% 
    dplyr::select(Parameter,Reproduction,Survival, `Reproduction + Survival`)  
  
  fit_table<-fit_table %>% 
    filter(Parameter %in% paramlatex$Parameter.latex)
  fit_table
}


table_wrangler_smart<-function(tables,col_headers,sigfigs,n_fixef, paramlatex){
  library(Matrix)
  library(tidyverse)

  output <- list()
  
  for(i in 1:length(tables)){
  table<-tables[[i]]

  
  table<-table %>% 
    mutate(across(where(is.numeric), ~ round(.x, digits = sigfigs))) %>% 
    mutate(Estimate = paste0(Estimate, " (", .lower,", ", .upper,")")) %>% 
    dplyr::select(Parameter, Estimate)
  
  fit_table<-table %>% 
    left_join(paramlatex) %>% 
    dplyr::select(-Parameter) %>% 
    rename(Parameter = Parameter.latex)%>% 
    dplyr::select(Parameter,Estimate)  
  
  fit_table<-fit_table %>% 
    filter(Parameter %in% paramlatex$Parameter.latex)
  output[[i]] <- fit_table
  
  colnames(output[[i]])[2] <- col_headers[i]
  
  }
  output
  output<-Reduce(full_join, output)
  output
  
  
  corr_order<-paramlatex %>% 
    filter(Parameter.latex %in% output$Parameter) %>% 
    dplyr::select(Parameter.latex) %>% 
    unique()

  
  output <-output %>% 
    arrange(factor(Parameter, levels = corr_order$Parameter.latex))
  
  newrow= c('\\textbf{Fixed effects} ($\\beta$)',NA, NA,NA)
  output <- rbind(newrow,output)
  
  
  r = n_fixef + 1 #one extra row since we've added the fixed effects header
  newrow= c('\\textbf{Random effects} ($\\sigma^2$)',NA, NA,NA)
  output<-rbind(output[1:r,],newrow,output[-(1:r),])
  
  output<-output %>% 
    rename('\\textbf{Parameter}' = Parameter)
  output
}



res_tables<-function(results_frame, table_capt, tbl_label,table_footnote, table_name,n_fixef){
  library(tables)
  library(knitr)
  library(kableExtra)
  
  options(knitr.kable.NA = '') #change NA values to blank

  kbl(results_frame,booktabs = TRUE,
      format = "latex",escape = FALSE, 
      caption = paste(table_capt), label = paste(tbl_label),digits = 2) %>% #
    #kable_styling(font_size = 12) %>%
    kable_classic(full_width = F, html_font = "Times New Roman") %>% 
    kable_styling(font_size = 10, position = "center") %>% 
    #kable_styling(latex_options="scale_down") %>% #fit to document margin
    pack_rows("Fixed effects $\\beta$", 1, n_fixef, escape = F) %>%
    pack_rows("Random effects $\\sigma^2$", n_fixef+1, nrow(results_frame), escape = F) %>% 
    add_footnote(table_footnote, notation = "symbol",escape = F) %>% 
    save_kable(paste0("Tables empirical data/Tables output/",table_name,".tex"))
}






stan_postplot<-function(stan_model,r0_fig_name,gamma_fig_name){
  
  library(bayesplot)
  
  library(ggplot2)
  
  
  elements<-row.names(summary(stan_model)$summary)
  
  isl_r0<-elements[which(str_detect(elements, "isl_r0"))]
  
  
  mcmc_areas(stan_model,point_est = "median", pars =c(print(isl_r0), "r_0"))+
    scale_y_discrete(
      labels = c("r_0" = "Meta-population",
                 "isl_r0[1]" = "Aldra",
                 "isl_r0[2]" = "Gjerøy",
                 "isl_r0[3]" = "Hestmannøy",
                 "isl_r0[4]" = "Indre Kvarøy",
                 "isl_r0[5]" = "Lovund",
                 "isl_r0[6]" = "Myken",
                 "isl_r0[7]" = "Nesøy",
                 "isl_r0[8]" = "Selvær",
                 "isl_r0[9]" = "Sleneset",
                 "isl_r0[10]" = "Træna"
      )
    )+
    xlab(bquote(r[0]))
  
  ggsave(paste0("Figures empirical data/",r0_fig_name), bg="white",width = 1950, height = 1365, units = "px")
  
  
  isl_g<-elements[which(str_detect(elements, "isl_g"))]
  
  
  mcmc_areas(stan_model,point_est = "median", pars =c(print(isl_g), "gamma"))+
    scale_y_discrete(
      labels = c("gamma" = "Meta-population",
                 "isl_g[1]" = "Aldra",
                 "isl_g[2]" = "Gjerøy",
                 "isl_g[3]" = "Hestmannøy",
                 "isl_g[4]" = "Indre Kvarøy",
                 "isl_g[5]" = "Lovund",
                 "isl_g[6]" = "Myken",
                 "isl_g[7]" = "Nesøy",
                 "isl_g[8]" = "Selvær",
                 "isl_g[9]" = "Sleneset",
                 "isl_g[10]" = "Træna"
                 
      )
    )+
    xlab(bquote(gamma))
  ggsave(paste0("Figures empirical data/",gamma_fig_name), bg="white", width = 1950, height = 1365, units = "px")
  
}





brm_postplot<-function(brm_model,density_variable,r0_fig_name,gamma_fig_name){
  library(tidybayes)
  
  library(ggplot2)

  
  if(density_variable == "n"){
    brm_islands<-brm_model %>% 
      spread_draws(b_Intercept, b_n,r_Location[Location,Parameter]) %>% 
      mutate(Estimate = case_when(Parameter == "n" ~ b_n+ r_Location,
                                  Parameter == "Intercept"~b_Intercept+r_Location))
  } else if(density_variable == "zwsc"){
    brm_islands<-brm_model %>% 
      spread_draws(b_Intercept, b_wsc_zpop,r_Location[Location,Parameter]) %>% 
      mutate(Estimate = case_when(Parameter == "wsc_zpop" ~ b_wsc_zpop + r_Location,
                                  Parameter == "Intercept"~b_Intercept+r_Location))
  } else {
    print("error, density_variable must be either n or zwsc")
  }


  r0_plot<-brm_islands %>%   
    filter(Parameter == "Intercept")
  
  
  r0_plot %>%   
    ggplot(aes(y = Location, x = Estimate)) +
    stat_halfeye(.width = c(0.5, 0.9), point_interval="median_qi",aes(fill = after_stat(level)))+
    scale_fill_brewer(na.translate = FALSE) +
    labs(y="",x=expression(r[0]))+
    theme_bw()+
    theme(legend.position = "none")+
    coord_cartesian(xlim = c(quantile(r0_plot$Estimate,probs =0.01),quantile(r0_plot$Estimate, probs = 0.99)))
  ggsave(paste0("Figures empirical data/",r0_fig_name), bg="white",width = 1950, height = 1365, units = "px")
  
  gamma_plot<-brm_islands %>%   
    filter(Parameter == "n"|Parameter == "wsc_zpop")
  
  gamma_plot %>%   
    ggplot(aes(y = Location, x = Estimate)) +
    stat_halfeye(.width = c(0.5, 0.9), point_interval="median_qi",aes(fill = after_stat(level)))+
    scale_fill_brewer(na.translate = FALSE) +
    labs(y="",x=expression(gamma))+
    theme_bw()+
    theme(legend.position = "none")+
    coord_cartesian(xlim = c(quantile(gamma_plot$Estimate,probs =0.01),quantile(gamma_plot$Estimate, probs = 0.99)))
  ggsave(paste0("Figures empirical data/",gamma_fig_name), bg="white",width = 1950, height = 1365, units = "px")
  
}

## Island predicted slopes at experienced pop size.  ####

brm_blups_islands<-function(brm_model,density_variable, fitness, divide_intercepts, average_metapars){
  library(tidyverse)
  library(tidybayes)


  outer<-brm_model$data %>% 
    filter(in_out == 1) %>% 
    dplyr::select(Location) %>% 
    unique()

  if(density_variable == "n"){

    brm_islands<-brm_model %>% 
      spread_draws(b_Intercept, b_n, b_in_out,`b_n:in_out`, r_Location[Location,Parameter]) %>% 
      mutate(Estimate = case_when(Parameter == "n" & !Location %in%  outer$Location ~ b_n+ r_Location,
                                  Parameter == "n" & Location %in%  outer$Location ~ b_n+ r_Location + `b_n:in_out`,
                                  Parameter == "Intercept" & !Location %in%  outer$Location ~ b_Intercept+r_Location,
                                  Parameter == "Intercept" & Location %in%  outer$Location ~b_Intercept+r_Location + b_in_out)) %>% 
      rename(Island = Location)
    
    if(average_metapars == "yes"){
      brm_islands<-brm_islands %>% 
        mutate(b_n = b_n + `b_n:in_out`/2) %>% 
        mutate(b_Intercept = b_Intercept + b_in_out/2)
    } else if(average_metapars == "no"){
      print("Warning, when metapars is set to no the metapopulation intercept and slope represent inner islands")
      brm_islands
    } else {
      print("error, average_metapars must be either yes or no")
    }
    
    if(divide_intercepts == "yes"){
      brm_islands<-brm_islands %>% 
        mutate(b_Intercept = b_Intercept - log(2)) %>% 
        mutate(Estimate = ifelse(Parameter == "Intercept", Estimate- log(2),
                                 Estimate))
    } else if(divide_intercepts == "no"){
      brm_full
    }  else {
      print("error, divide_intercepts must be either yes or no")
    }
    
  
    
    brm_metapop<-brm_islands %>% 
      ungroup() %>% 
      rename(r0 = b_Intercept, g = b_n) %>% 
      median_qi(r0, g) %>%
      mutate(Island = "meta-population") %>% 
      dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
    
    brm_blups<- brm_islands %>% 
      median_qi(Estimate) %>% 
      pivot_wider(values_from = c(Estimate, .lower, .upper), names_from = Parameter) %>% 
      rename(r0 = Estimate_Intercept , g = Estimate_n , r0.lower = .lower_Intercept, r0.upper = .upper_Intercept, g.lower = .lower_n, g.upper = .upper_n) %>% 
      dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
    
    }else if(density_variable == "N"){
      
      brm_islands<-brm_model %>% 
        spread_draws(b_Intercept, b_N,b_in_out,`b_N:in_out`, r_Location[Location,Parameter]) %>% 
        mutate(Estimate = case_when(Parameter == "N" & !Location %in%  outer$Location ~ b_N+ r_Location,
                                    Parameter == "N" & Location %in%  outer$Location ~ b_N+ r_Location + `b_N:in_out`,
                                    Parameter == "Intercept" & !Location %in%  outer$Location ~ b_Intercept+r_Location,
                                    Parameter == "Intercept" & Location %in%  outer$Location ~b_Intercept+r_Location + b_in_out)) %>% 
        rename(Island = Location)
      
      
      if(average_metapars == "yes"){
        brm_islands<-brm_islands %>% 
          mutate(b_N = b_N + `b_N:in_out`/2) %>% 
          mutate(b_Intercept = b_Intercept + b_in_out/2)
      } else if(average_metapars == "no"){
        print("Warning, when metapars is set to no the metapopulation intercept and slope represent inner islands")
        brm_islands
      } else {
        print("error, average_metapars must be either yes or no")
      }
      
      
      if(divide_intercepts == "yes"){
        brm_islands<-brm_islands %>% 
          mutate(b_Intercept = b_Intercept - log(2)) %>% 
          mutate(Estimate = ifelse(Parameter == "Intercept", Estimate- log(2),
                                   Estimate))
      } else if(divide_intercepts == "no"){
        brm_full
      }  else {
        print("error, divide_intercepts must be either yes or no")
      }
      

      
      brm_metapop<-brm_islands %>% 
        ungroup() %>% 
        rename(r0 = b_Intercept, g = b_N) %>% 
        median_qi(r0, g) %>%
        mutate(Island = "meta-population") %>% 
        dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
      
      brm_blups<- brm_islands %>% 
        median_qi(Estimate) %>% 
        pivot_wider(values_from = c(Estimate, .lower, .upper), names_from = Parameter) %>% 
        rename(r0 = Estimate_Intercept , g = Estimate_N , r0.lower = .lower_Intercept, r0.upper = .upper_Intercept, g.lower = .lower_N, g.upper = .upper_N) %>% 
        dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
      
    }else if(density_variable == "mc_N"){
      
      brm_islands<-brm_model %>% 
        spread_draws(b_Intercept, b_mc_N, b_in_out,`b_mc_N:in_out`, r_Location[Location,Parameter]) %>% 
        mutate(Estimate = case_when(Parameter == "mc_N" & !Location %in%  outer$Location ~ b_mc_N+ r_Location,
                                    Parameter == "mc_N" & Location %in%  outer$Location ~ b_mc_N+ r_Location + `b_mc_N:in_out`,
                                    Parameter == "Intercept" & !Location %in%  outer$Location ~ b_Intercept+r_Location,
                                    Parameter == "Intercept" & Location %in%  outer$Location ~b_Intercept+r_Location + b_in_out)) %>% 
        rename(Island = Location)
      
      if(average_metapars == "yes"){
        brm_islands<-brm_islands %>% 
          mutate(b_mc_N = b_mc_N + `b_mc_N:in_out`/2) %>% 
          mutate(b_Intercept = b_Intercept + b_in_out/2)
      } else if(average_metapars == "no"){
        print("Warning, when metapars is set to no the metapopulation intercept and slope represent inner islands")
        brm_islands
      } else {
        print("error, average_metapars must be either yes or no")
      }
      
      
      if(divide_intercepts == "yes"){
        brm_islands<-brm_islands %>% 
          mutate(b_Intercept = b_Intercept - log(2)) %>% 
          mutate(Estimate = ifelse(Parameter == "Intercept", Estimate- log(2),
                                   Estimate))
      } else if(divide_intercepts == "no"){
        brm_full
      }  else {
        print("error, divide_intercepts must be either yes or no")
      }
      
 
      
      brm_metapop<-brm_islands %>% 
        ungroup() %>% 
        rename(r0 = b_Intercept, g = b_mc_N) %>% 
        median_qi(r0, g) %>%
        mutate(Island = "meta-population") %>% 
        dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
      
      brm_blups<- brm_islands %>% 
        median_qi(Estimate) %>% 
        pivot_wider(values_from = c(Estimate, .lower, .upper), names_from = Parameter) %>% 
        rename(r0 = Estimate_Intercept , g = Estimate_mc_N , r0.lower = .lower_Intercept, r0.upper = .upper_Intercept, g.lower = .lower_mc_N, g.upper = .upper_mc_N) %>% 
        dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
    
  } else if(density_variable == "zmc_N"){

    brm_islands<-brm_model %>% 
      spread_draws(b_Intercept, b_zNR,b_in_out,`b_zNR:in_out`, r_Location[Location,Parameter]) %>% 
      mutate(Estimate = case_when(Parameter == "zNR" & !Location %in%  outer$Location ~ b_zNR+ r_Location,
                                  Parameter == "zNR" & Location %in%  outer$Location ~ b_zNR+ r_Location + `b_zNR:in_out`,
                                  Parameter == "Intercept" & !Location %in%  outer$Location ~ b_Intercept+r_Location,
                                  Parameter == "Intercept" & Location %in%  outer$Location ~b_Intercept+r_Location + b_in_out))%>% 
      rename(Island = Location)
    
    
    if(average_metapars == "yes"){
      brm_islands<-brm_islands %>% 
        mutate(b_zNR = b_zNR + `b_zNR:in_out`/2) %>% 
        mutate(b_Intercept = b_Intercept + b_in_out/2)
    } else if(average_metapars == "no"){
      print("Warning, when metapars is set to no the metapopulation intercept and slope represent inner islands")
      brm_islands
    } else {
      print("error, average_metapars must be either yes or no")
    }
    
    if(divide_intercepts == "yes"){
      brm_islands<-brm_islands %>% 
        mutate(b_Intercept = b_Intercept - log(2)) %>% 
        mutate(Estimate = ifelse(Parameter == "Intercept", Estimate- log(2),
                                 Estimate))
    } else if(divide_intercepts == "no"){
      brm_full
    }  else {
      print("error, divide_intercepts must be either yes or no")
    }
    

    
    brm_metapop<-brm_islands %>% 
      ungroup() %>% 
      rename(r0 = b_Intercept, g = b_zNR) %>% 
      median_qi(r0, g) %>%
      mutate(Island = "meta-population") %>% 
      dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
    
    brm_blups<- brm_islands %>% 
      median_qi(Estimate) %>% 
      pivot_wider(values_from = c(Estimate, .lower, .upper), names_from = Parameter) %>% 
      rename(r0 = Estimate_Intercept , g = Estimate_zNR , r0.lower = .lower_Intercept, r0.upper = .upper_Intercept, g.lower = .lower_zNR, g.upper = .upper_zNR) %>% 
      dplyr::select(Island, r0, r0.lower, r0.upper, g, g.lower, g.upper)
  } else {
    print("error, density_variable must be either n, N, mc_N or zmc_N")
  }
  

  
  brm_full<- brm_blups %>% 
    bind_rows(brm_metapop)
  
  brm_full<-brm_full %>% 
    mutate_if(is.character, str_replace_all, pattern = fixed("."), replacement = " ") #fixed() is needed around the dot otherwise the dot is interpreted as "all text".
  brm_full
}
##Extract fixed parameter estimates for islands modelled as separate populations. results_single_islands is a list of model outputs ####

get_single_island_estimates<-function(results_single_islands,divide_intercepts,summary){

  ind_islands<-names(results_single_islands)
  
  all_isl_res <- NULL;
  
  for(i in 1:length(ind_islands)){
    
    single_island <- results_single_islands[[i]]
    
    single_isl_draws<-single_island %>% 
      spread_draws(b_Intercept, b_gamma) %>% 
      rename(r0 = b_Intercept) %>% 
      rename(g = b_gamma)
    
    if(divide_intercepts == "yes"){
      single_isl_draws <- single_isl_draws %>% 
        mutate(r0 = r0-log(2))
    } else if(divide_intercepts == "no"){
      single_isl_draws
    }  else {
      print("error, divide_intercepts must be either yes or no")
    }
    
    if(summary == FALSE){
    single_isl_post_dist<-single_isl_draws %>% 
      dplyr::select(r0, g)
    single_isl_post_dist$Island <- names(results_single_islands)[i]
    all_isl_res <- bind_rows(all_isl_res, single_isl_post_dist)
    
    } else if(summary == TRUE){
    single_isl_draws<-single_isl_draws %>%
      median_qi(r0,g)
    
    single_isl_draws$Island <- names(results_single_islands)[i]
    single_isl_draws<-single_isl_draws %>% 
      relocate(Island, .before = r0) %>% 
      dplyr::select(-.width, -.point, -.interval)
    
    all_isl_res <- bind_rows(all_isl_res, single_isl_draws)
    }
    else{
      print("ERROR! summary must be TRUE or FALSE")
    }
  }
  all_isl_res
  

}


#Prob defines the interval of experienced pop size for the v on n plot. A prob of
#0.8 exludes the lower and upper 10% quantiles of the experienced pop sizes for the
#predicted values of individual islands.

island_slopes<-function(post_frame,popsize_frame,data_observed,mean_centred, density_variable, prob,v_ymin, v_ymax,w_ymin, w_ymax, fitness){
  library(tidyverse)
  library(rcartocolor)
  n_isl <- post_frame %>% 
    filter(Island != "meta-population") %>% 
    nrow()
  
  safe_colourblind_palette <- carto_pal(n_isl, "Safe")
  
  if("meta-population" %in% post_frame$Island){
  post_frame<-post_frame %>% 
    mutate(metapop = ifelse(Island == "meta-population", 1,
                            0)) %>% 
    group_by(metapop) %>% 
    arrange(metapop,Island, by_group = T) %>% 
    ungroup()
  } else {
    post_frame<-post_frame %>% 
    arrange(Island)
  }
  

  
  #If the predictor is mean centred, we start and stop the x-axis at the deviation furthest from 0, rather than from 0 to max popsize.

  if(mean_centred == "yes"){
    v_xlim = seq(min(popsize_frame$popsize), max(popsize_frame$popsize), length.out = 1000)
    } else if(mean_centred == "no"){
    v_xlim = seq(0,max(popsize_frame$popsize), length.out=1000)
  
  }else{
    print("Error, mean_centred must be yes or no!")
  }
  
  


  par(mfrow=c(1,2))
  par(mar = c(10,10,10,10))
  par(mgp=c(1,2,0)) #change distance betweens axis tick lables and tick marks. Firsta value is y, second x, and the third is the line holding the tick marks
  global_par<-post_frame %>% 
    filter(Island == "meta-population")
  
  island_par<-post_frame %>% 
    filter(Island != "meta-population")
  if("meta-population" %in% post_frame$Island){
  island_colours<-c(safe_colourblind_palette, "#000000")
  } else{
    island_colours<-safe_colourblind_palette
  }
  
  plot(1, type = "n", ann = FALSE,
       cex.axis = 3,
       xlim = c(min(v_xlim), max(v_xlim)), 
       ylim = c(v_ymin, v_ymax))
  
  if(density_variable == "n"){
    x_label_v <- "Log adult density n"
  } else if(density_variable == "mc_N"){
    x_label_v <- expression("N"[mu])
  } else if(density_variable == "zmc_N"){
    xLab <- "Mean-centred and scaled adult density"
    x_label_v <- bquote(.(xLab[1])~"N"["m"])
   } else if(density_variable == "N"){
      x_label_v <- expression("N")
  } else{
    print("error, density_variable must be either n,N,zmc_N or mc_N")
  }

  title(xlab = x_label_v, cex.lab = 3.5,
        line = 8)
  if(fitness == "r2s"){
    y_label_v <- "log average individual demographic contribution"
    #y_label_v <- bquote(.(y_label_v[1])~"v"["rS"])
  }else if(fitness == "recruitment"){
    y_label_v = "log recruit production per female"
    #y_label_v <- bquote(.(y_label_v[1])~"v"["r"])
    
  }else if(fitness == "survival"){
    y_label_v <- expression("v"["S"])
  } else{
    print("error, fitness must be either r2s,recruitment or survival")
  }
 
  title(ylab = y_label_v, cex.lab = 3.5, line = 6)


  
  for(i in 1:length(island_par$Island)){
    
    ind_island<-subset(popsize_frame, Island == island_par$Island[i])
    popsize_isl <- seq(quantile(ind_island$popsize,((1-prob)/2)),quantile(ind_island$popsize,1-((1-prob)/2)), length.out=1000)
    ind_island_par = subset(island_par, Island == island_par$Island[i])
  
    v_isl<-ind_island_par$r0 + ind_island_par$g*popsize_isl
   
 
    
    points(v_isl~popsize_isl, col=island_colours[i], lty = 1, type ="l", lwd = 7)
   
    }
    
  

  if("meta-population" %in% post_frame$Island){
  v_glo<-global_par$r0 + global_par$g*v_xlim
  points(v_glo ~ v_xlim, col = "black", type = "l", lwd = 11)
  } else {}
  
  plot(1, type = "n", ann = FALSE,
       cex.axis = 3,
       xlim = c(0, max(data_observed$estimated_pop_size)), 
       ylim = c(w_ymin,w_ymax))
  

  title(xlab = "Adult density N", cex.lab = 3.5,
        line = 8)
  
  if(fitness == "r2s"){
    y_label_w <- "Average individual demographic contribution"
    #y_label_w <- bquote(.(y_label_w[1])~"w"["rS"])
  }else if(fitness == "recruitment"){
    y_label_w <- "Recruit production per female"
    #y_label_w <- bquote(.(y_label_w[1])~"w"["r"])
    
  }else if(fitness == "survival"){
    y_label_w <- expression("w"["S"])
  } else{
    print("error, fitness must be either r2s,recruitment or survival")
  }
  
  
  title(ylab = y_label_w, cex.lab = 3.5,
        line = 6)
  
  
  if(density_variable == "n"){
    N_metapop = exp(v_xlim)
  }else{
  N_metapop <- seq(min(data_observed$estimated_pop_size), max(data_observed$estimated_pop_size), length.out = 1000)
  }


  for(i in 1:length(island_par$Island)){
    
    ind_island<-subset(popsize_frame, Island == island_par$Island[i])
    popsize_isl <- seq(min(ind_island$popsize),max(ind_island$popsize), length.out=1000)
    ind_island_par = subset(island_par, Island == island_par$Island[i])
    obs_pop_isl<- subset(data_observed,Location == island_par$Island[i])
    N_popsize_isl <- seq(min(obs_pop_isl$estimated_pop_size), max(obs_pop_isl$estimated_pop_size), length.out = 1000)
    
    v_isl<-ind_island_par$r0 + ind_island_par$g*popsize_isl
    
    if(density_variable == "n"){
      plot_popsize_isl = exp(popsize_isl)
    } else{
      plot_popsize_isl = N_popsize_isl
    }

    
 
    points(exp(v_isl)~plot_popsize_isl, col=island_colours[i], lty = 1, type ="l", lwd = 7)
    
  }
  if("meta-population" %in% post_frame$Island){
  points(exp(v_glo) ~ N_metapop, col = "black", type = "l", lwd = 11)
  } else{}
  
  if("meta-population" %in% post_frame$Island){
    meta_lty = 1
    
  } else{
    meta_lty = NA
  }
  
  legend("topright", legend = str_to_title(post_frame$Island),
         
         col=c(island_colours), lty=1,lwd=5, cex=2.5,
         title="Island", text.font=4, bg='darkgrey') #lty=c(1:11,meta_lty) if we want different line types for all the lines.

  
}

island_slopes_single<-function(post_frame, popsize_frame, mean_centred, prob,w_ymin, w_ymax, fitness){
  library(tidyverse)
  library(rcartocolor)
  n_isl <- nrow(post_frame) %>% 
    filter(Island != "meta-population") %>% 
    
  
  safe_colourblind_palette <- carto_pal(n_isl, "Safe")
  
  
  if("meta-population" %in% post_frame$Island){
    post_frame<-post_frame %>% 
      mutate(metapop = ifelse(Island == "meta-population", 1,
                              0)) %>% 
      group_by(metapop) %>% 
      arrange(metapop,Island, by_group = T) %>% 
      ungroup()
  } else {
    post_frame<-post_frame %>% 
      arrange(Island)
  }
  
  #If the predictor is mean centred, we start and stop the x-axis at the deviation furthest from 0, rather than from 0 to max popsize
  if(mean_centred == "yes"){
    plot_popsizes = seq(min(popsize_frame$popsize),max(popsize_frame$popsize), length.out=1000)
    
  } else{plot_popsizes = seq(0,max(popsize_frame$popsize), length.out=1000)
  }
  
  
  par(mfrow=c(1,1))
  par(mar = c(10,10,10,10))
  par(mgp=c(1,2,0)) #change distance betweens axis tick lables and tick marks. Firsta value is y, second x, and the third is the line holding the tick marks
  global_par<-post_frame %>% 
    filter(Island == "meta-population")
  
  island_par<-post_frame %>% 
    filter(Island != "meta-population")
  
  if("meta-population" %in% post_frame$Island){
    island_colours<-c(safe_colourblind_palette, "#000000")
  } else{
    island_colours<-safe_colourblind_palette
  }
  
  plot(1, type = "n", ann = FALSE,
       cex.axis = 3,
       xlim = c(min(plot_popsizes), max(plot_popsizes)), 
       ylim = c(w_ymin, w_ymax))
  
  
  title(xlab = "n", cex.lab = 3.5,
        line = 8)
  if(fitness == "r2s"){
    y_label_v <- expression("w"["rS"])
  }else if(fitness == "recruitment"){
    y_label_v <- expression("w"["r"])
    
  }else if(fitness == "survival"){
    y_label_v <- expression("w"["S"])
  } else{
    print("error, fitness must be either r2s,recruitment or survival")
  }
  
  title(ylab = y_label_v, cex.lab = 3.5, line = 6)
  
  
  
  
  for(i in 1:length(island_par$Island)){
    
    ind_island<-subset(popsize_frame, Location == island_par$Island[i])
    popsize_isl <- seq(quantile(ind_island$popsize,((1-prob)/2)),quantile(ind_island$popsize,1-((1-prob)/2)), length.out=1000)
    
    
    ind_island_par = subset(island_par, Island == island_par$Island[i])
    
    
    v_isl<-ind_island_par$r0 + ind_island_par$g*popsize_isl
    
    
    
    points(exp(v_isl)~popsize_isl, col=island_colours[i], lty = 1, type ="l", lwd = 7)
    
  }
  legend("topright", legend = str_to_title(post_frame$Island),
         
         col=c(island_colours), lty=1,lwd=5, cex=2.5,
         title="Island", text.font=4, bg='darkgrey')#lty=c(1:11) if we want different line types for the legend
}


