data {
  // Number of clusters (an integer)
  int<lower=0> nind;            // Ind ID and number of individuals, since there is only one individual per row.
  int<lower=2> n_occasions;     // Number of capture occasions (years)
  int<lower=0> nflok;
  int<lower=0> nFY;
    // Clusters indentifiers
  int<lower=1,upper=nflok> flok[nind]; //Islands. Only changes among individuals (no dispersal as adults)
  int<lower=1,upper=nFY> FY[nflok, n_occasions-1];    //island-years
  //Predictors
  matrix[nflok, n_occasions-1] n; //density
  matrix[nind, n_occasions - 1] age; //This is mean centred age, so it is not an integer. changes among individuals and among years.
  matrix[nind, n_occasions - 1] age_sq; // //slope, quadratic age, added 6/11-2023 at 09:00
  int<lower=0,upper=1> sex[nind];  //changes among individuals
  int<lower=0, upper = 1> in_out[nind]; // inner or outer system, changes among islands but we specify is as the length of individuals, just like flok.
  //Capture history
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
}
//
transformed data {
  // Compoud declaration is enabled in Stan 2.13
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];
  int<lower=0,upper=n_occasions> last[nind];
//
  //  n_occ_minus_1 = n_occasions - 1;
  for (i in 1:nind) {
    for (k in 1:n_occasions)
      if (y[i,k]) {
        first[i] = k;
        break;
      }
    for (k_rev in 1:n_occasions) {
      int k = n_occasions - k_rev + 1;
      if (y[i,k]) {
        last[i] = k;
        break;
      }
    }
  }
}
//
parameters {
  real<lower=0,upper=1> mean_phi;    // Mean survival
  real<lower=0,upper=1> mean_p;      // Mean recapture
  real gamma; //slope, density reg
   real B_age; //slope, age
   real B_age_sq;//slope, quadratic age
  real B_sex; //slope, sex
  real B_in_out;//effects of inner vs outer system slope
  real interac_in_out_gamma; //interaction gamma and inner outer
  real interac_age_sex;
  //sds
  vector<lower=0>[1] sigmaFY_p;  //sd island-year recapture prob
  vector<lower=0>[1] sigmaI_phi; //sd individual mean survival
  vector<lower=0>[1] sigmaFY_phi; //sd year-island mean survival
  vector<lower=0>[2] sigmaF_phi; //sd island mean survival, intercept and slope
  //random effects
  matrix[nind,1] zI; //vector scaled individual blups (intercepts). These have to be of type matrix because you can't multiply two vectors in stan
  matrix[nFY,1] zFY; 
  matrix[2,nflok] zF; //matrix scaled island blups(intercepts and slopes)
  cholesky_factor_corr[2] LF;  // factor to estimate covariance int-slopes for islands (Flok). To estimate covariance. Here we specify the off-diagonal number of cells in the variance covariance matrix. We have two random effect levels, so there are two.
  matrix[nFY, 1] zp; //matrix scaled island year recapture probability
}
//
transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi; //survival
  matrix<lower=0,upper=1>[nflok, n_occ_minus_1] p; //recapture prob
  //
  vector[nind] epsilonI_phi; //random effect of individual survival
  vector[nFY] epsilonFY_phi; //random effect of individual survival
  matrix[2,nflok] epsilonF_phi; //random effect of island survival
  vector[nFY] epsilonFY_p; //island year recpture prob
//
epsilonI_phi  =  zI * sigmaI_phi; // ind intercept
epsilonFY_phi  =  zFY * sigmaFY_phi; // year-island intercept
epsilonF_phi = diag_pre_multiply(sigmaF_phi, LF) * zF;// equation for random Location intercept and slope
epsilonFY_p = zp * sigmaFY_p; //year-island recapture prob
//    
//
  real mu_p; //logit transformed mean recapture prob
  real mu_phi; //logit transformed mean survival
  //
//
  mu_p = logit(mean_p);
  mu_phi = logit(mean_phi);
  // Constraints
  for (i in 1:nind) {
    phi[i, 1:(first[i] - 1)] = rep_row_vector(0,first[i] - 1);
    for (t in first[i]:n_occ_minus_1) {
      phi[i, t] = inv_logit(mu_phi+ epsilonI_phi[i] + epsilonF_phi[1,flok[i]] +epsilonFY_phi[FY[flok[i],t]] + (gamma+epsilonF_phi[2,flok[i]])*n[flok[i],t] + B_sex*sex[i] + B_age*age[i, t]+ B_age_sq * age_sq[i,t] + B_in_out*in_out[i] + interac_in_out_gamma*in_out[i]*n[flok[i],t] +interac_age_sex*age[i,t]*sex[i]);
    }
  }
  for (j in 1:nflok)
    for (t in 1:n_occ_minus_1)
      p[j, t] = inv_logit(mu_p + epsilonFY_p[FY[j, t]]);
}
model {
  // Priors
  // Uniform priors are implicitly defined.
  mean_phi ~ beta(2, 2); 
  mean_p ~ beta(10, 3); 
  gamma ~ normal(0,2); 
  B_age ~ normal(0,2);
  B_age_sq ~ normal(0,2);
  B_sex ~ normal(0,2);
  B_in_out ~ normal(0,2);
  interac_in_out_gamma ~ normal(0,2);
  interac_age_sex ~ normal(0,2);
  to_vector(zI) ~ normal(0, 1); //this is scaling the sampling of the blups to have mean 0 and sd 1. This is a z standardisation
  to_vector(sigmaI_phi) ~ normal(0, 1); 
  to_vector(zF) ~ normal(0, 1);
  to_vector(sigmaF_phi) ~ normal(0, 1); 
  to_vector(zFY) ~ normal(0, 1);
  to_vector(sigmaFY_phi) ~ normal(0, 1); 
  to_vector(zp) ~ normal(0,1);
  to_vector(sigmaFY_p) ~ normal(0, 1); 
  // Likelihood
  for (i in 1:nind) {                     //loop over individuals to assess survival probability
    int j = flok[i];                      //this varies by island and year. The specific island an individual is in is specified here
    for (t in first[i]:(last[i] - 1)) {   //here, we say that between the first capture occasion and the last capture occasion the probability of an individual being alive is 1 (but shouldn't the loop go from first + 1 since it ends at last - 1?)
      1 ~ bernoulli(phi[i, t]);          //here we say that it is alive between first and last time it was seen
      y[i, t + 1] ~ bernoulli(p[j, t]); //needs to be + 1 because survival is transition between years and recapture rates are defined by the year following the transition.
    }
    if(last[i]<n_occasions) {                  //after the last time an individual was seen it is much trickier to know if an individual was alive and not recapture or if it was dead. Therefore, we need to assess the probability that an individual was alive but not captured. This we can only do for years before the last time step, because after the last time step we don't have any more information.
      real chi = 1;                            //chi is the probability that an individual was alive but not seen. We set a starting value of chi as 1. This does not mean that the probability of not being captured is 1, it is just a value to stick in the equation below. Since it is 1 it doesn't change the results of the equation for the first round. However, for the second round of the loop, the new chi calculated from the first round is used. If phi is 0.5 and p is 0.8, chi of the first round is 0.6. When this is used in the equation in the second round, the probability will be multiplied by 0.6. The new chi will thus be even smaller. This new, smaller chi is then used to multiply the equation at the third round. Therefore the probability of being alive but not captured is reduced by each year. Do this for long enough then the probability of being alive but not captured is so small that the individual is probably dead. 
      for (t in 1:(n_occasions-last[i])) {     //stan cannot loop from the back to the front, so we need to do a trick. By starting a loop at t (which starts at 1) to the number of occassions (which is not the actual calendar year it is the total number of years) minus the last time an individual was seen
        int t_curr = n_occasions - t;           // we can create a vector of current time and therefore a loop that goes from the last occasion of the time series to the last time an individual was seen.
        chi = (1 - phi[i,t_curr]) + phi[i,t_curr] * (1 - p[j,t_curr]) * chi; //we calculate the probability of not being observed as the sum of the probability of being dead (1- phi) and the probaility of being alive but not captured (phi * (1-p)). For the first round, this is multiplied to chi = 1. But a new chi is generated here, which is used to multiply the equation for the next round. So, each round (year) chi is getting smaller
      }
      1 ~ bernoulli(chi);
    }
  }
}
generated quantities{
real<lower=0> Sigma2_phi; // Individual intercept variance. Transform sd to var
//real<lower=0> Sigma2_gamma; // Individual slope variance. Transform sd to var
real<lower=0> Sigma2_Isl_phi; // Island intercept variance. Transform sd to var
real<lower=0> Sigma2_Isl_gamma; // Island slope variance. Transform sd to var
real<lower=0> Sigma2_YI_phi; // Year_island slope variance. Transform sd to var
real<lower=0> Sigma2_YI_p; // Year_island recapture probability
real cov_Isl; // the covariance between r0 and gamma, islands
row_vector[nflok] isl_phi; //island r0s. Must be row_vector because we specified Isl as a matrix with BLUP's on the rows.
row_vector[nflok] isl_g; //island gammas. Must be row_vector because we specified Isl as a matrix with BLUP's on the rows.
//
//
//matrix[2, 2]  Omega_I;
matrix[2, 2]  Omega_Isl;
Sigma2_phi=sigmaI_phi[1]^2; //this gives the variance of the individual level intercepts
//
Sigma2_Isl_phi= sigmaF_phi[1]^2; // this gives the variance of the island level intercepts
Sigma2_Isl_gamma= sigmaF_phi[2]^2; // this gives the variance of the island level slopes
Sigma2_YI_phi= sigmaFY_phi[1]^2; //this gives the variance of the Year_island intercepts
Sigma2_YI_p = sigmaFY_p[1]^2; //island-year recapture prob variance
isl_phi=mu_phi+epsilonF_phi[1,]; //add global r0 to island r0 BLUP's.
isl_g= gamma + epsilonF_phi[2,]; //add global gamma to island gamma BLUP's
//
Omega_Isl = LF * LF'; //this gives the correlation matrix for islands. Multiply L by the inverse of L (L')
cov_Isl = Omega_Isl[1,2]*sqrt(Sigma2_Isl_gamma*Sigma2_Isl_phi);// backtransform the correlation (off diagonal of Omega matrix) to covariance for the island level
}
