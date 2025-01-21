#A possible function to allow the option for multiple point estimate methods!


brm2table_short<-function(brm_result,type, width){
  library(tidybayes)
  library(brms)
  library(tidyverse)

  
  post_est<-function(x,type, width){
    if(type == "median")(
      return(median_qi(x, .width = width)))
    if(type == "mean")(
      return(mean_qi(x, .width = width)))
    if(type == "mode")(
      return(mode_qi(x, .width = width)))
    else(
      print("Error! type must be mean, median or mode!"))
    
  }

  if(brm_result$family[1] == "negbinomial")(
  draws<-brm_result %>% 
    gather_draws(`b_.*`,shape, regex = TRUE)
  
  )else (  draws<-brm_result %>% 
             gather_draws(`b_.*`, regex = TRUE))
  
    brm_result_res <- post_est(draws,type, width) %>% 
    rename(Parameter = .variable, Estimate=.value)
  brm_result_res
  
  #The variance component point estimates from VarCorr are the mean. We want median.
  # we can do this by extracting the posterior distribution and do median_qi instead.
  #To get a correct variance estimate for components only available as sd we need to get
  #the posterior distribution, square all values, and then calculate the point estimate.
  #Just taking the square of a point estimate changes the outcome alot because of rounding errors.
  varcovs <- list()



  varcovs[[1]]<-post_est(as.data.frame((VarCorr(brm_result,summary = FALSE)$ID$sd)^2), type, width) #just taking the square of the sd estimate can create so many rounding errors that the estimate is unreliable. Here we take the square of the posterior distribution of sd instead.  
  varcovs[[2]]<-post_est(as.data.frame(VarCorr(brm_result,summary = FALSE)$Location$cov)[,1],type, width)
  varcovs[[3]]<-post_est(as.data.frame(VarCorr(brm_result,summary = FALSE)$Location$cov)[,4], type, width)
  varcovs[[4]] <-post_est(as.data.frame((VarCorr(brm_result, summary = FALSE)$Year_Location$sd)^2),type, width) #just taking the square of the sd estimate can create so many rounding errors that the estimate is unreliable. Here we take the square of the posterior distribution of sd instead.  
  varcovs[[5]]= post_est(as.data.frame(VarCorr(brm_result,summary = FALSE)$Location$cov)[,2],type, width)
  
  
  
  varnames <- c("Sigma2_r0","Sigma2_Isl_r0", "Sigma2_Isl_gamma", "Sigma2_YI", "cov_Isl")
  
  varcovs_table<-NULL;
  for(i in 1:5){
  colnames(varcovs[[i]])[1] <- "Estimate"
  colnames(varcovs[[i]])[2] <- ".lower"
  colnames(varcovs[[i]])[3] <- ".upper"
  rownames(varcovs[[i]]) <- NULL
  varcovs[[i]]$Parameter <- varnames[i]
  varcovs_table<- rbind.data.frame(varcovs_table, varcovs[[i]])

  }
  varcovs_table
 
  
 
  
  brm_result_res<-brm_result_res %>% 
    bind_rows(varcovs_table)
  brm_result_res
}
