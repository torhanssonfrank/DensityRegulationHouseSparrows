data {
  // Number of clusters (an integer)
  int<lower=0> nind;            // Ind ID and number of individuals, since there is only one individual per row.
  int<lower=0> n_occasions;     // Number of capture occasions (years)
  int<lower=0> nflok;
  int<lower=0> nFY;
  int<lower=0> max_startyear;
    // Clusters indentifiers
  int<lower=1,upper=nflok> flok[nind,n_occasions-1]; //Islands. changes within individuals
  int<lower=1,upper=nFY> FY[nflok, n_occasions-1];    //island-years
  int<lower=1,upper=max_startyear> start_year[nflok]; //the year when the data series starts for each island
  //real<lower=0>  max_age;                                  // Maximum age
  //Predictors
  matrix[nflok, n_occasions-1] obs_pop; //observed population sizes
  matrix[nflok, n_occasions-1] obs_rec; //observed number of recruits. flok is the born island, not obs island
  matrix[nind, n_occasions-1] NR; //nestling to recruit survival design matrix
  matrix[nind, n_occasions-1] AD; //Adult to adult survival design matrix
  matrix[nind, n_occasions - 1] age; //This is mean centred age, so it is not an integer. changes among individuals and among years. before specified as int<lower=-1, upper=max_age> age[nind, n_occasions - 1]
  matrix[nind, n_occasions - 1] age_sq; // mean centred age squared
  int<lower=0,upper=1> sex[nind];  //changes among individuals
  int<lower=0, upper = 1> in_out[nflok]; // inner or outer system, changes among islands but we specify is as the length of individuals, just like flok.
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
  real<lower=0,upper=1> mean_phi_rec;      // Mean survival nestlings to rec
  real<lower=0,upper=1> mean_p;      // Mean recapture
  real B_in_out;
  real B_in_out_rec;
  real gamma;
  real gamma_rec;
  real gamma_in_out;
  real gamma_in_out_rec;
  real B_sex;
  real B_age;
  real B_age_sq;
  real B_sex_age;
  //
   vector<lower=0>[1] sigmaFY_p;  //sd island-year recapture prob
   matrix[nFY, 1] zp; //matrix scaled island survival probability rec
   //FY phi
   vector<lower=0>[1] sigmaFY_phi;  //sd island survival prob adults
   vector<lower=0>[1] sigmaFY_phi_rec;  //sd island survival prob recruits
   matrix[nFY, 1] zFYphi; //matrix scaled island year survival probability ads
   matrix[nFY, 1] zFYphi_rec; //matrix scaled island year survival probability rec
   //F phi
   vector<lower=0>[2] sigmaF_phi;
   matrix[2, nflok] zphi; //matrix scaled island survival probability ads
   cholesky_factor_corr[2] LF; 
   //
   vector<lower=0>[2] sigmaF_phi_rec;
   matrix[2, nflok] zphi_rec; //matrix scaled island survival probability rec
   cholesky_factor_corr[2] LF_rec; 
}
//
transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi; //survival
  matrix<lower=0,upper=1>[nflok, n_occ_minus_1] p; //n_occ_minus_1 recapture prob
  matrix[nflok, n_occasions-1] N_corr;
  matrix[nflok, n_occasions-1] Rec_corr;
  matrix[nflok, n_occasions-1] zNR;
//
vector[nFY] epsilonFY_p; //island year recpture prob
//
vector[nFY] epsilonFY_phi; //island year survival prob ads
vector[nFY] epsilonFY_phi_rec; //island year survival prob recs
//
matrix[2,nflok] epsilonF_phi; //random effect of island survival
matrix[2,nflok] epsilonF_phi_rec;
//
epsilonFY_p = zp * sigmaFY_p; //year-island recapture prob
//
epsilonFY_phi = zFYphi * sigmaFY_phi; //
epsilonFY_phi_rec = zFYphi_rec * sigmaFY_phi_rec; //
//
epsilonF_phi = diag_pre_multiply(sigmaF_phi, LF) * zphi;// equation for random Location intercept and slope
epsilonF_phi_rec =diag_pre_multiply(sigmaF_phi_rec, LF_rec) * zphi_rec; //
//
  real mu_p; //logit transformed mean recapture prob
  real mu_phi; //logit transformed mean survival
  real mu_phi_rec;
  //
//
  mu_p = logit(mean_p);
  mu_phi = logit(mean_phi);
  mu_phi_rec = logit(mean_phi_rec);
  // Constraints
    for (j in 1:nflok)
    for (t in 1:n_occ_minus_1)
      p[j, t] = inv_logit(mu_p + epsilonFY_p[FY[j, t]]);  //
for (j in 1:nflok){
  N_corr[j,1:(start_year[j]-1)] = rep_row_vector(0,start_year[j]-1); //add 0 for islands that don't have data at the start of the study period
  Rec_corr[j,1:(start_year[j]-1)] = rep_row_vector(0,start_year[j]-1);
  //
  N_corr[j, start_year[j]] = obs_pop[j, start_year[j]]/mean_p; //start_year is the first year when we have nestling data for an island.
  Rec_corr[j, start_year[j]] = obs_rec[j, start_year[j]]/mean_p;
  for (t in (start_year[j]+1):n_occ_minus_1){
  N_corr[j, t] = obs_pop[j, t] / p[j, t-1];
  Rec_corr[j, t] = obs_rec[j, t] / p[j, t-1];
  }
  zNR[j,1:(start_year[j]-1)] = rep_row_vector(0,start_year[j]-1);
  zNR[j,start_year[j]:n_occ_minus_1] = (N_corr[j,start_year[j]:n_occ_minus_1] - mean(N_corr[j,start_year[j]:n_occ_minus_1]))./40;
}
  for (i in 1:nind) {
    phi[i, 1:(first[i] - 1)] = rep_row_vector(0,first[i] - 1);
    for (t in first[i]:n_occ_minus_1) {
      int j = flok[i,t]; //individuals can change island
      //when we include sex it is super important that we have an adult design matrix! Otherwise all juveniles will be females
      phi[i, t] = inv_logit(mu_phi + AD[i,t]* (epsilonF_phi[1,j] + epsilonFY_phi[FY[j, t]] +  B_in_out*in_out[j]
      + B_sex * sex[i] + B_age * age[i, t] + B_age_sq * age_sq[i, t] + B_sex_age * sex[i] * age[i, t] ) //phi adults
      //
      + NR[i,t]* (mu_phi_rec + epsilonF_phi_rec[1,j] + epsilonFY_phi_rec[FY[j, t]] + B_in_out_rec*in_out[j])  //phi recruits
      //
      + zNR[j, t]* NR[i,t]*(gamma_rec + epsilonF_phi_rec[2,j] + gamma_in_out_rec * in_out[j]) // gamma recruits
      //
      + zNR[j, t]* AD[i, t]* (gamma+ epsilonF_phi[2,j] +  gamma_in_out * in_out[j])); //gamma adults
      //
    }
  }
}
model {
  // Priors
  mean_phi ~ uniform(0, 1); //uniform
  mean_phi_rec ~  uniform(0, 1); //uniform
  mean_p ~ beta(1.5, 1.2); //beta prior. weakly informative
  B_in_out ~ normal(0, 2);
  B_in_out_rec ~ normal(0, 2);
  B_sex ~ normal(0, 2);
  B_age ~ normal(0, 2);
  B_age_sq ~ normal(0, 2);
  B_sex_age ~ normal(0, 2);
  gamma ~ normal(0, 2);
  gamma_rec ~ normal(0, 2); 
  gamma_in_out ~ normal(0, 2);
  gamma_in_out_rec ~ normal(0, 2);
  to_vector(zp) ~ normal(0,1);
  to_vector(sigmaFY_p) ~ normal(0, 1); //
  //
  to_vector(zFYphi) ~ normal(0,1);
  to_vector(sigmaFY_phi) ~ normal(0, 1);
  to_vector(zFYphi_rec) ~ normal(0,1);
  to_vector(sigmaFY_phi_rec) ~ normal(0, 1);
  //
  to_vector(zphi) ~ normal(0,1);
  to_vector(sigmaF_phi) ~ normal(0, 1);
  LF ~ lkj_corr_cholesky(4); //  weakly informative
  LF_rec ~ lkj_corr_cholesky(4); //  weakly informative
  //
  to_vector(zphi_rec) ~ normal(0,1);
  to_vector(sigmaF_phi_rec) ~ normal(0, 1);
  // Likelihood
  for (i in 1:nind) {  
    for (t in first[i]:(last[i] - 1)) {   //here, we say that between the first capture occasion and the last capture occasion the probability of an individual being alive is 1. Survival is estimated forwards in time. Therefore, the last capture occasion assesses survival from the second last capture occasion to the last.
      1 ~ bernoulli(phi[i, t]);          //here we say that it is alive between first and last time it was seen
      int j = flok[i,t];                      //this varies by island and year. The specific island an individual is in is specified here
      y[i, t + 1] ~ bernoulli(p[j,t]); //Here we assess recapture rates between the first and last time an individual was seen. THis needs to be + 1 because survival is transition between years and recapture rates are defined by the year following the transition.
    }
    if(last[i]<n_occasions) {                  //after the last time an individual was seen it is much trickier to know if an individual was alive and not recapture or if it was dead. Therefore, we need to assess the probability that an individual was alive but not captured. This we can only do for years before the last time step, because after the last time step we don't have any more information.
      real chi = 1;                            //chi is the probability that an individual was not seen. We set a starting value of chi as 1. This does not mean that the probability of not being captured is 1, it is just a value to stick in the equation below. Since it is 1 it doesn't change the results of the equation for the first round. However, for the second round of the loop, the new chi calculated from the first round is used. If phi is 0.5 and p is 0.8, chi of the first round is 0.6. When this is used in the equation in the second round, the probability will be multiplied by 0.6. The new chi will thus be even smaller. This new, smaller chi is then used to multiply the equation at the third round. Therefore the probability of being alive but not captured is reduced by each year. Do this for long enough then the probability of being alive but not captured is so small that the individual is probably dead. 
       int l = flok[i,last[i]];                 //after the last time an individual was seen we only look in the last island it was seen.
      for (t in 1:(n_occasions-last[i])) {     //stan cannot loop from the back to the front, so we need to do a trick. By starting a loop at t (which starts at 1) to the number of occassions (which is not the actual calendar year it is the total number of years) minus the last time an individual was seen
        int t_curr = n_occasions - t;           // we can create a vector of current time and therefore a loop that goes from the last occasion of the time series to the last time an individual was seen.
        chi = (1 - phi[i,t_curr]) + phi[i,t_curr] * (1 - p[l,t_curr]) * chi; //we calculate the probability of not being observed as the sum of the probability of being dead (1- phi) and the probaility of being alive but not captured (phi * (1-p)). For the first round, this is multiplied to chi = 1. But a new chi is generated here, which is used to multiply the equation for the next round. So, each round (year) chi is getting smaller
      }
      1 ~ bernoulli(chi);
    }
  }
}
generated quantities{
  //estimate number of nestlings
    matrix[nflok, n_occasions-1] Nest_corr;
    vector[nflok] Nest_corr_F;
    //
      for (j in 1:nflok){
    Nest_corr_F[j] = sum(Rec_corr[j, (start_year[j]+1):n_occ_minus_1])/inv_logit(mu_phi + mu_phi_rec + B_in_out_rec*in_out[j] + epsilonF_phi_rec[1,j]); //Adult blups not needed when we have adult design matrix
    //
    Nest_corr[j,1:(start_year[j]-1)] = rep_row_vector(0,start_year[j]-1); //add 0 for islands that don't have data at the start of the study period
    Nest_corr[j, start_year[j]] = Rec_corr[j, start_year[j]]/inv_logit(mu_phi + mu_phi_rec + B_in_out_rec*in_out[j] + epsilonF_phi_rec[1,j]);
    for (t in (start_year[j]+1):n_occ_minus_1){
    //
  Nest_corr[j, t] = Rec_corr[j, t]/inv_logit(mu_phi + mu_phi_rec + B_in_out_rec*in_out[j] + epsilonF_phi_rec[1,j] 
                                              + epsilonFY_phi_rec[FY[j, t-1]]
                                              + (gamma_rec + gamma_in_out_rec*in_out[j]
                                              + epsilonF_phi_rec[2,j])*zNR[j, t-1]); //
}
}
//
//Create projection matrices and calculate elasticities
  vector<lower=0>[nflok] f;
  vector<lower=0>[nflok] f_simple;
  vector[nflok] s_n;
  vector[nflok] s_a;
for (j in 1:nflok){
  f[j] =  (sum(Nest_corr[j, (start_year[j]+1):(n_occasions-1)])/2)/(sum(N_corr[j, start_year[j]:(n_occasions-2)])/2);//divide by 2 to get female recruits and female nestlings. Since we do it to both the denominator and the numerator it doesn't make a difference. It is just to be explicit about female demographic dominance.
  f_simple[j] =  (Nest_corr_F[j]/2)/(sum(N_corr[j, start_year[j]:(n_occasions-2)])/2);////divide by 2 to get female recruits and female nestlings
  s_n[j] = inv_logit(mu_phi+  mu_phi_rec+B_in_out_rec*in_out[j] + epsilonF_phi_rec[1,j]);
  s_a[j] = inv_logit(mu_phi+ B_in_out*in_out[j] +epsilonF_phi[1,j]);
}
matrix[2, 2] A[nflok];
  real lambda[nflok];
  vector[2] w[nflok]; // right eigenvectors (stable stage)
  vector[2] v[nflok]; // left eigenvectors (reproductive value)
  matrix[2, 2] elasticity[nflok];
//
  real tol = 1e-8;
//
  for (l in 1:nflok) {
    // 1. Build projection matrix A for population l
    A[l,1,1] = 0;
    A[l,1,2] = f_simple[l];      // nestling production. 
    A[l,2,1] = s_n[l];    // nestling-to-adult survival
    A[l,2,2] = s_a[l];    // adult survival
//
    // 2. Power iteration for right eigenvector (w[l])
    vector[2] w_old = rep_vector(1.0, 2);
    for (i in 1:100) {
      vector[2] w_new = A[l] * w_old;
      w_new /= sum(w_new); // normalize
      if (max(abs(w_new - w_old)) < tol) break;
      w_old = w_new;
    }
    w[l] = w_old;
//
    // 3. Dominant eigenvalue ??[l] (finite rate of increase)
    vector[2] Aw = A[l] * w[l];
    lambda[l] = dot_product(w[l], Aw) / dot_product(w[l], w[l]);
//
    // 4. Power iteration for left eigenvector (v[l])
    vector[2] v_old = rep_vector(1.0, 2);
    for (i in 1:100) {
      vector[2] v_new = A[l]' * v_old;
      v_new /= sum(v_new);
      if (max(abs(v_new - v_old)) < tol) break;
      v_old = v_new;
    }
    v[l] = v_old;
//
    // 5. Normalize left eigenvector so that dot(v, w) = 1
    real norm_factor = dot_product(v[l], w[l]);
    v[l] /= norm_factor;
//
    // 6. Compute elasticity matrix
    for (i in 1:2)
      for (j in 1:2)
        elasticity[l,i,j] = (A[l,i,j] / lambda[l]) * v[l][i] * w[l][j];
  }
  //random intercepts
  real<lower=0> Sigma2_YI_p;
  real<lower=0> Sigma2_FY_phi;
  real<lower=0> Sigma2_FY_phi_rec;
  real<lower=0> Sigma2_isl_phi;
  real<lower=0> Sigma2_isl_phi_rec;
  //random slopes
  real<lower=0> Sigma2_isl_gamma;
  real<lower=0> Sigma2_isl_gamma_rec;
   //covariance
  real cov_isl; //
  matrix[2, 2]  Omega_isl;
  //
  real cov_isl_rec; //
  matrix[2, 2]  Omega_isl_rec;
  //
  Sigma2_YI_p = sigmaFY_p[1]^2;
  Sigma2_FY_phi = sigmaFY_phi[1]^2;
  Sigma2_FY_phi_rec = sigmaFY_phi_rec[1]^2;
  //
  Sigma2_isl_phi = sigmaF_phi[1]^2;
  Sigma2_isl_gamma = sigmaF_phi[2]^2;
  //
  Sigma2_isl_phi_rec = sigmaF_phi_rec[1]^2;
  Sigma2_isl_gamma_rec = sigmaF_phi_rec[2]^2;
  //
  Omega_isl = LF * LF';
  cov_isl = Omega_isl[1,2]*sqrt(Sigma2_isl_gamma*Sigma2_isl_phi); //Omega_isl[1,2] holds the correlation
  //
  Omega_isl_rec = LF_rec * LF_rec';
  cov_isl_rec = Omega_isl_rec[1,2]*sqrt(Sigma2_isl_gamma_rec*Sigma2_isl_phi_rec); 
}
