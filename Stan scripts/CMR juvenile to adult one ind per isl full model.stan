data {
  // Number of clusters (an integer)
  int<lower=0> nind;            // Ind ID and number of individuals, since there is only one individual per row.
  int<lower=0> n_occasions;     // Number of capture occasions (stages)
  int<lower=0> nflok;
  int<lower=0> nFY;
    // Clusters identifiers
  int<lower=0, upper = nflok> flok[nind];
  int<lower=0, upper=nFY> FY[nind];
  //Predictors
  vector[nind] N; //density
  vector[nind] in_out; //inner outer system
  matrix[nind, n_occasions] Fl;
  matrix[nind, n_occasions] R; //design matrix recruits.
  matrix[nind, n_occasions] A; //design matrix adults
  // design matrices recapture rates
  matrix[nind, n_occasions] Fl_p; //design matrix fledglings recap
  matrix[nind, n_occasions] R_p; //design matrix recruits recap
  matrix[nind, n_occasions] A_p; //design matrix recruits recap
  //Capture history
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
}
//
transformed data {
  // Compoud declaration is enabled in Stan 2.13
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> last[nind];
// Find the last occasion a focal individual was seen. First is the same for all
  for (i in 1:nind) {
    for (k_rev in 1:n_occasions) { //trick to loop backwards
      int k = n_occasions - k_rev + 1;
      if (y[i,k]) { //y[i,k] is the same as y[i,k] == 1 in Stan
        last[i] = k;
        break;
      }
    }
  }
}
//
parameters {
  real<lower=0,upper=1> mean_phi;    // Mean survival (fledgelings)
  real<lower=0,upper=1> mean_rec_phi;    // Mean survival recruits
  real<lower=0,upper=1> mean_ad_phi;    // Mean survival adults
  real<lower=0,upper=1> mean_p;      // Mean recapture (fledgelings)
  real<lower=0,upper=1> mean_rec_p;      // Mean recapture (recruits)
  real<lower=0,upper=1> mean_ad_p;      //Mean recapture adults
  real gamma;
  real gamma_rec;
  real phi_fledge_interac_inout;
  real phi_rec_interac_inout;
  real gamma_fledge_interac_inout;
  real gamma_rec_interac_inout;
  //sds intercepts and slopes survival density reg
  vector<lower=0>[1] sigmaFY_phi_fledge;//stage specific year island random intercepts fledge
  vector<lower=0>[1] sigmaFY_phi_rec;//stage specific year island random intercepts rec
  vector<lower=0>[2] sigmaF_fledge; //sd island mean survival, intercepts and slopes
  vector<lower=0>[2] sigmaF_rec; //sd island mean survival, intercepts and slopes
  //recapture sds
  vector<lower=0>[1] sigmaFY_p_fledge;//stage specific year island random intercepts for recapture probability ledge
  vector<lower=0>[1] sigmaFY_p_rec;//stage specific year island random intercepts for recapture probability recruits
  //random effects survival and density reg 
  matrix[nFY,1] zFY_fledge;
  matrix[nFY,1] zFY_rec;
  matrix[2,nflok] zFledge;//matrix scaled island blups(intercepts and slopes) fledge
  matrix[2,nflok] zRec;
  cholesky_factor_corr[2] L_fledge; // factor to estimate covariance int-slopes fledglings
  cholesky_factor_corr[2] L_rec; // factor to estimate covariance int-slopes recruits
  //random effects recapture prob
  matrix[nFY,1] zFY_p_fledge;
  matrix[nFY,1] zFY_p_rec;
}
//
transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi; //survival 
  matrix<lower=0,upper=1>[nind, n_occasions] p; //recapture probability
  //year island blups
  vector[nFY] epsilonFY_phi_fledge; //random intercept year island survival fledge
  vector[nFY] epsilonFY_phi_rec; //random intercept year island survival rec
  //island blups
  matrix[2,nflok] epsilonF_fledge; //random intercepts and slopes fledglings island
  matrix[2,nflok] epsilonF_rec; //random intercepts and slopes recruits island
  //recapture blups
 vector[nFY] epsilonFY_p_fledge; //random intercept year island recapture fledge
 vector[nFY] epsilonFY_p_rec; //random intercept year island recapture recruits
//
epsilonFY_phi_fledge =  zFY_fledge * sigmaFY_phi_fledge; // year-island intercept fledge
epsilonFY_phi_rec =  zFY_rec * sigmaFY_phi_rec; // year-island intercept fledge
//
epsilonF_fledge = diag_pre_multiply(sigmaF_fledge, L_fledge)* zFledge;
//
epsilonF_rec = diag_pre_multiply(sigmaF_rec, L_rec)* zRec;
//
epsilonFY_p_fledge = zFY_p_fledge * sigmaFY_p_fledge;
epsilonFY_p_rec = zFY_p_rec * sigmaFY_p_rec;
//
  real mu_p; //logit transformed mean recapture prob
  real mu_rec_p; //logit transformed recruit recapture prob
  real mu_ad_p; //logit transformed adult recapture prob
  real mu_phi; //logit transformed mean survival
  real mu_rec_phi; //
  real mu_ad_phi; //
  //
//
  mu_p = logit(mean_p);
  mu_rec_p = logit(mean_rec_p);
  mu_ad_p = logit(mean_ad_p);
  mu_phi = logit(mean_phi);
  mu_rec_phi = logit(mean_rec_phi);
  mu_ad_phi = logit(mean_ad_phi);
  // Constraints
  for (i in 1:nind) {
    for (t in 1:n_occ_minus_1) {
    phi[i, t] = inv_logit(mu_phi + (epsilonF_fledge[1,flok[i]] + 
    epsilonFY_phi_fledge[FY[i]])*Fl[i,t]+ 
    //
    phi_fledge_interac_inout*in_out[i]*Fl[i,t] +
    //
    (mu_rec_phi + epsilonF_rec[1,flok[i]] + epsilonFY_phi_rec[FY[i]]) * R[i,t] + 
    //
    phi_rec_interac_inout*in_out[i]*R[i,t] +
    //
    mu_ad_phi * A[i,t] +
    //
    (gamma + epsilonF_fledge[2,flok[i]])*Fl[i,t]*N[i]+
    //
    gamma_fledge_interac_inout*in_out[i]*Fl[i,t]*N[i] +
    //
    (gamma_rec + epsilonF_rec[2,flok[i]]) * R[i,t] *N[i] +
    //
    gamma_rec_interac_inout*in_out[i]*R[i,t]*N[i]);//
    //
    }
     p[i,1] = 1;
    for (t in 2:n_occasions) { 
        p[i,t] = inv_logit(mu_p +
    (epsilonFY_p_fledge[FY[i]]) * Fl_p[i,t] + 
     (mu_rec_p + epsilonFY_p_rec[FY[i]]) * R_p[i,t] +
     mu_ad_p*A_p[i,t]); 
    }
  }
     }
model {
  // Priors
  // Uniform priors are implicitly defined.
  mean_phi ~ normal(0, 2); //
  mean_rec_phi ~ normal(0, 2); //
  mean_ad_phi ~ normal(0, 2); //
  phi_fledge_interac_inout ~ normal(0, 2);
  phi_rec_interac_inout ~ normal(0, 2);
  //
  gamma ~ normal(0, 2);
  gamma_rec ~ normal(0, 2);
  gamma_fledge_interac_inout ~ normal(0, 2);
  gamma_rec_interac_inout ~ normal(0, 2);
  //
  mean_p ~ beta(10, 3); //I assume that fledgling recapture is the same as recruits
  mean_rec_p ~ beta(10, 3);
  mean_ad_p ~ uniform(0, 1); //there is a lot of heterogeneity here so we specify a uniform prior
  //
  to_vector(epsilonFY_p_fledge) ~ normal(0,1); 
  to_vector(epsilonFY_p_rec) ~ normal(0,1); 
  to_vector(sigmaFY_phi_fledge) ~ normal(0,1);
  to_vector(sigmaFY_phi_rec) ~ normal(0,1);
  to_vector(sigmaF_fledge) ~ normal(0,1);
  to_vector(sigmaF_rec) ~ normal(0,1);
  to_vector(zFY_fledge) ~ normal(0,1);
  to_vector(zFY_rec) ~ normal(0,1);
  to_vector(zFledge) ~ normal(0,1);
  to_vector(zRec) ~ normal(0,1);
  to_vector(zFY_p_fledge) ~ normal(0,1);
  to_vector(zFY_p_rec) ~ normal(0,1);
  L_fledge ~ lkj_corr_cholesky(4); // unclear why there is a 4.
  L_rec ~ lkj_corr_cholesky(4); // unclear why there is a 4.
  // Likelihood
   for (i in 1:nind) {
    int j = FY[i]; 
    for (t in 2:(last[i])) {
      1 ~ bernoulli(phi[i, t-1]);
      y[i, t] ~ bernoulli(p[j, t]);
    }
    if(last[i]<n_occasions) {
      real chi = 1;
      for (t in 1:(n_occasions-last[i])) {
        int t_curr = n_occasions - t;
        int t_next = t_curr + 1;
        chi = (1 - phi[i,t_curr]) + phi[i,t_curr] * (1 - p[j,t_next]) * chi;
      }
      1 ~ bernoulli(chi);
    }
  }
}
generated quantities{
//variance components
//
real<lower=0> Sigma2_YI_phi_fledge; //year island survival fledgelings
real<lower=0> Sigma2_YI_phi_rec; //year island survival recruits
//
real<lower=0> Sigma2_isl_phi_fledge; // island survival fledge
real<lower=0> Sigma2_isl_phi_rec; // island survival rec
//
real<lower=0> Sigma2_isl_g_fledge; // island density reg fledge
real<lower=0> Sigma2_isl_g_rec; // island density reg rec
//
real cov_fledge;
real cov_rec;
matrix[2, 2]  Omega_fledge;
matrix[2, 2]  Omega_rec;
//
real<lower=0> Sigma2_YI_p_fledge; //year island recapture prob fledge
real<lower=0> Sigma2_YI_p_rec; //year island recapture prob recruits
//BLUP's
row_vector[nflok] isl_phi_fledge;
row_vector[nflok] isl_phi_rec;
//
row_vector[nflok] isl_g_fledge;
row_vector[nflok] isl_g_rec;
//
Sigma2_YI_phi_fledge = sigmaFY_phi_fledge[1]^2;
Sigma2_YI_phi_rec = sigmaFY_phi_rec[1]^2;
//
Sigma2_isl_phi_fledge = sigmaF_fledge[1]^2;
Sigma2_isl_phi_rec = sigmaF_rec[1]^2;
//
Sigma2_isl_g_fledge = sigmaF_fledge[2]^2;
Sigma2_isl_g_rec = sigmaF_rec[2]^2;
//
Omega_fledge = L_fledge * L_fledge'; //this gives the correlation matrix for islands. Multiply L by the inverse of L (L')
cov_fledge = Omega_fledge[1,2]*sqrt(Sigma2_isl_g_fledge*Sigma2_isl_phi_fledge);// backtransform the correlation (off diagonal of Omega matrix) to covariance for the island level
Omega_rec = L_rec * L_rec';
cov_rec = Omega_rec[1,2]*sqrt(Sigma2_isl_g_rec*Sigma2_isl_phi_rec);
//
Sigma2_YI_p_fledge = sigmaFY_p_fledge[1]^2;
Sigma2_YI_p_rec = sigmaFY_p_rec[1]^2;
//
isl_phi_fledge=mu_phi+epsilonF_fledge[1,]; //add intercept to island phi BLUP's
isl_phi_rec=mu_phi+mu_rec_phi+epsilonF_rec[1,]; //add intercept to rec phi beta and to island phi BLUP's
//
isl_g_fledge=gamma+epsilonF_fledge[2,]; //add intercept to island phi BLUP's
isl_g_rec=gamma+gamma_rec+epsilonF_rec[2,]; //add intercept to rec phi beta and to island phi BLUP's
}
