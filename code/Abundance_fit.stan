// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/aseed_sets/case_studies/bayeseed_sparse_regreseed_sion.html

data{
  int<lower = 1> N; // Number of observations within each year
  int<lower = 1> S; // Number of plant species (same acroseed_s years for consistency)
  
  matrix[N,S] SpAbundance;  // Fecundity of the focal species in each observation
  real<lower=0>  growth_ratio[N];  // Fecundity of the focal species in each observation
  real<lower=0> lambda_mean_obs;
  
  int<lower = 1> year[N]; // Indicator variable for the year each observations
  int<lower = 1> Y; // number of years
  

  real<lower=0> g_mean; // germination
  real<lower=0> g_sd; // germination
  real<lower=0> seed_s_mean; // seed survival
  real<lower=0> seed_s_sd; // seed survival

  matrix[1,S] alpha_initial_mean;
  matrix[1,S] alpha_slope_mean;
  matrix[1,S] c_mean;
  matrix[1,S] N_opt_mean;
  matrix[1,S] alpha_initial_sd;
  matrix[1,S] alpha_slope_sd;
  matrix[1,S] c_sd;
  matrix[1,S]  N_opt_sd;
  
}

parameters{
  real<lower=0> lambda_mean[1];
  real<lower=0> lambda_sd[Y];
  real<lower=0> g[1];
  real<lower=0> seed_s[1];
  
  real<lower=-1,upper=1> alpha_initial[S];

  real<lower=-1,upper=0> alpha_slope[S];

  real<lower=-1,upper=0> c[S];

  real<lower = 0> N_opt[S];
  
  vector<lower=0>[1] disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{

  // Declare objects neceseed_sary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N] interaction_effects;
  
  // loop parameters
  vector[N] lambda_ei;
  matrix[N,S] alpha_value;
  

 // implement the biological model
  for(n in 1:N){
    lambda_ei[n] = lambda_mean[1] + lambda_sd[year[n]];
    for(s in 1:S){
    alpha_value[n,s] = alpha_initial[s] + (c[s] * (1 - exp( alpha_slope[s] * (SpAbundance[n,s]-N_opt[s]))))/(1+exp(alpha_slope[s] * (SpAbundance[n,s]-N_opt[s])));
    }
    
    interaction_effects[n] = sum(alpha_value[n,] .* SpAbundance[n,]);
    
    }
}
model{
  for(s in 1:S){
  N_opt[s] ~ normal(1, N_opt_sd[1,s]);
  alpha_initial[s] ~ normal(alpha_initial_mean[1,s], alpha_initial_sd[1,s]);  // normal prior on lambda
  alpha_slope[s] ~ normal(alpha_slope_mean[1,s], alpha_slope_sd[1,s]);     // normal prior on alpha
  c[s] ~ normal(c_mean[1,s], c_sd[1,s]);     // normal prior on alpha
    }
  disp_dev ~ normal(0,1);
  g ~ normal(g_mean, g_sd);  
  seed_s ~ normal(seed_s_mean, seed_s_sd); 
  lambda_mean ~ normal(0,1);
  for(y in 1:Y){
    lambda_sd[y] ~ normal(0,1);
  }
  
  for(n in 1:N){
  growth_ratio[n] ~ normal((1 - g[1])*seed_s[1] + g[1]*lambda_ei[n]*exp(interaction_effects[n]),disp_dev);
  }
}
