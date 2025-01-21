#Goodnes of fit ####





#Sensitivities ####

growth_factors<-growth_factors %>% 
  filter(lambda<4)

growth_factors<-growth_factors %>% 
  mutate(n = log(N))

hist(growth_factors$lambda, breaks = 100)
descdist(growth_factors$lambda, boot = 500, discrete = F)


N_mod <- glmer(lambda ~ N +(1|Year) + (N|Location), 
               family = Gamma(link = "log"),
               data = growth_factors)

n_mod <- glmer(lambda ~ n + (1|Year) + (n|Location),
               family = Gamma(link = "log"),
               data = growth_factors)

AIC(N_mod,n_mod)