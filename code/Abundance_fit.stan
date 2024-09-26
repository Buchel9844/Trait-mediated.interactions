// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/aseed_sets/case_studies/bayeseed_sparse_regreseed_sion.html

data{
  int<lower = 1> N; // Number of observation total
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int focal_pos; //position of intra in SpAbundance dataset 
  int<lower = 1> numb_year; // number of years
  int<lower = 1> numb_plot; // number of plot observation within each year
 
  matrix[N,S] SpAbundance;  // Abundance of the species in each observation 
  int obs_spI[N]; //vector of obs for focal
  int year[N];  // year vector for each obs
  int plot[N];  //plot vector for each obs

  real<lower=0> g_mean; // germination
  real<lower=0> g_sd; // germination
  real<lower=0> seed_s_mean; // seed survival
  real<lower=0> seed_s_sd; // seed survival
  
  real PDSI_mean[numb_year];
  real PDSI_sd[numb_year];

  matrix[1,S] alpha_initial_mean;
  matrix[1,S] alpha_slope_mean;
  matrix[1,S] c_mean;
  matrix[1,S] N_opt_mean;
  matrix[1,S] alpha_initial_sd;
  matrix[1,S] alpha_slope_sd;
  matrix[1,S] c_sd;
  matrix[1,S] N_opt_sd;
  
}

parameters{
  
  real<lower=0> g[1];
  real<lower=0> seed_s[1];
  real PDSI[numb_year];

  real<lower=-1,upper=1> alpha_initial[S];

  real<lower=-1,upper=0> alpha_slope[S];

  real<lower=-1,upper=0> c[S];

  real<lower = 0> N_opt[S];
  
  vector<lower=0>[1] disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  real beta[1];
  real<lower = 0> lambda_init[1];
}

transformed parameters{

  // Declare objects neceseed_sary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N] interaction_effects;
  
  vector<lower = 0>[numb_year] lambda_ei; // loop parameter
  vector[N] seedtotal; // estimated number of seed per year and plot
  vector[S] alpha_value;
 // implement the biological model
    for(n in 1:N){

    lambda_ei[year[n]] = lambda_init[1] + beta[1]*PDSI_mean[year[n]];
    
    //print("plot",plot[n]);
    
    //print("year",year[n]);
    
    for(s in 1:S){
    alpha_value[s] = alpha_initial[s] + (c[s] * (1 - exp( alpha_slope[s] * (SpAbundance[n,s]-N_opt[s]))))/(1+exp(alpha_slope[s] * (SpAbundance[n,s]-N_opt[s])));
    }
   
   interaction_effects[n] = sum(alpha_value * SpAbundance[n,]); //Abundance_n[y,p] call the matrix of year y, and row p
  
    seedtotal[n] = (1-g[1])*seed_s[1]*obs_spI[n]./g[1] + obs_spI[n]*lambda_ei[year[n]]*exp(interaction_effects[n]);
    //print(seedtotal[n]);
    }
}

model{
  vector[N] stems_spI; // estimated number of stem per year and plot

  for(s in 1:S){
  N_opt[s] ~ normal(N_opt_mean[1,s], N_opt_sd[1,s]);
  alpha_initial[s] ~ normal(alpha_initial_mean[1,s], alpha_initial_sd[1,s]);  // normal prior on lambda
  alpha_slope[s] ~ normal(alpha_slope_mean[1,s], alpha_slope_sd[1,s]);     // normal prior on alpha
  c[s] ~ normal(c_mean[1,s], c_sd[1,s]);     // normal prior on alpha
    }
    
  disp_dev ~ normal(0,1);
  seed_s ~ normal(seed_s_mean, seed_s_sd); 
  beta  ~ normal(0,1);
  lambda_init ~ normal(0,1);
  g ~ normal(g_mean, g_sd);  

   for(n in 1:N){
  //PDSI[year[n]] ~ normal(PDSI_mean[year[n]],PDSI_sd[year[n]]); // draw a value of PDSI proper to that year

  stems_spI[n] = seedtotal[n]*g[1]; // get the number of stem for the year y and plot p for species I 
  // the probability distribtuion of species I for year y follows a poisson distribution include all observation across plots
  //print(stems_spI[n]);
  if(stems_spI[n] > 0)
  obs_spI[n] ~ poisson(stems_spI[n]); // draw a value of observed species I for that year 
  }
}

generated quantities{
  vector[N] stems_spI; // estimated number of stem per year and plot
  vector[N] pred_spI;
  for(n in 1:N){
    stems_spI[n] = seedtotal[n]*g[1];
  if(stems_spI[n] <= 0) break ;
    pred_spI[n] = poisson_lpmf(obs_spI[n]|stems_spI[n]); // draw a value of observed species I for that year 
  }
}
