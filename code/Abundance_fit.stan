// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/aseed_sets/case_studies/bayeseed_sparse_regreseed_sion.html

data{
  int<lower = 1> N; // Number of observation total
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int<lower = 1> numb_year; // number of years
  int<lower = 1> numb_plot; // number of plot observation within each year
  real ratio_scale;
  
  matrix[N,S] SpAbundance;  // Abundance of the species in each observation 
  matrix[numb_year,1] SpAbundance_med;
  int obs_spI[N]; //vector of obs for focal
  int year[N];  // year vector for each obs
  int plot[N];  //plot vector for each obs

  real PDSI_mean[numb_year];
  //real<lower=0> g_vec[numb_year]; // germination

  //real<lower=0> seed_s; // seed survival
  real lambda_mean_prior;
  real lambda_sd_prior;
  
  //real beta_prior;
  real beta_low;
  real beta_up;
  real beta_prior;
  real alpha_initial_mean[S];
  real alpha_slope_mean[S];
  real  c_mean[S];
  real N_opt_mean[S];
  
}

parameters{
  real<lower=beta_low,upper=beta_up> beta[1];
  real<lower=-1,upper=1> g_beta[1];
  real<lower=0,upper=1> g_init[1];
   real<lower=0,upper=1> seed_s[1];
  vector<lower=0>[1] lambda_mean;
}

transformed parameters{

  // Declare objects for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N-numb_plot] interaction_effects;
    matrix[N-numb_plot,S] alpha_value;
    matrix[N-numb_plot,S] alpha_ei;
    vector[N-numb_plot] seedtotal; // estimated number of seed per year and plot
    vector<lower = 0>[N-numb_plot] lambda_ei; // loop parameter
    vector<lower = 0>[N] g_ei; // loop parameter
     matrix[N-numb_plot,S] SpAbundance_scaled;
   for(n in 1:N){
     g_ei[n] = g_init[1] + g_beta[1]*exp(PDSI_mean[year[n]]);
   }
   for(n in 1:N-numb_plot){
    lambda_ei[n] = lambda_mean[1]+ beta[1]*exp(PDSI_mean[year[n]]);
  
   for(s in 1:S){
    SpAbundance_scaled[n,s] = SpAbundance[n,s]*ratio_scale;// 15/25
    alpha_value[n,s] = alpha_initial_mean[s] + (c_mean[s] * (1 - exp(alpha_slope_mean[s] * (SpAbundance_scaled[n,s]-N_opt_mean[s]))))/(1+exp(alpha_slope_mean[s] * (SpAbundance_scaled[n,s]- N_opt_mean[s])));
    alpha_ei[n,s]=alpha_value[n,s]*SpAbundance_scaled[n,s];
    }
    
   interaction_effects[n] = sum(alpha_ei[n,]);
  if (obs_spI[n] > 0){
    //g_vec[year[n]]
   seedtotal[n] = obs_spI[n] *((1- g_ei[n])*seed_s[1]./g_ei[n] + lambda_ei[n]*exp(interaction_effects[n]));
   }
  else{
   seedtotal[n] = SpAbundance_med[year[n],1]*((1- g_ei[n])*seed_s[1]./g_ei[n] +  lambda_ei[n]*exp(interaction_effects[n]));
   }
  }
}

model{
  matrix[1,N-numb_plot] stems_spI ; // estimated number of stem per year and plot
  seed_s ~  beta(20,20); //normal(seed_s_mean,seed_s_sd); //beta(5,5);
  beta  ~ normal(beta_prior,1); //normal(beta_prior,1)
  g_beta  ~ normal(0,0.1); //normal(beta_prior,1)
  g_init  ~  beta(10,10); //normal(beta_prior,1)
  lambda_mean ~ normal(lambda_mean_prior,lambda_sd_prior);
  
 for(n in 1:N-numb_plot){
  stems_spI[1,n] = seedtotal[n]*g_ei[n+numb_plot]; // get the number of stem for the year y and plot p for species I 
  
  obs_spI[n+numb_plot] ~  poisson(stems_spI[1,n]); // draw a value of observed species I for that year 

  }
}
