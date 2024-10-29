// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/aseed_sets/case_studies/bayeseed_sparse_regreseed_sion.html

data{
  int<lower = 1> N; // Number of observation total
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int<lower = 1> numb_year; // number of years
  int<lower = 1> numb_plot; // number of plot observation within each year

  matrix[N,S] SpAbundance;  // Abundance of the species in each observation 
  matrix[numb_year,1] SpAbundance_med;
  int obs_spI[N]; //vector of obs for focal
 // real obs_spI[N]; //vector of obs for focal
  int year[N];  // year vector for each obs
  int plot[N];  //plot vector for each obs

  real PDSI_mean[numb_year];
  real<lower=0> g_mean; // germination
  real<lower=0> g_sd_up; // germination
    real<lower=0> g_sd_min; // germination

  real<lower=0> seed_s_mean; // seed survival
  real<lower=0> seed_s_sd_up; // seed survival
    real<lower=0> seed_s_sd_min; // seed survival

  real beta_prior;
  real beta_prior_low;
  real beta_prior_up;
  real alpha_initial_mean[S];
  real alpha_slope_mean[S];
  real  c_mean[S];
  real N_opt_mean[S];
  real lambda_mean_prior;
   real lambda_sd_prior;
  
}

parameters{
  
  real<lower=seed_s_sd_min,upper=1> seed_s[1];

  real<lower=0> disp_dev[1]; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  //real<lower = beta_prior_low,upper=beta_prior_up> beta[1];
  real<lower = -20,upper=20> beta[1];
   // real<lower = 0> lambda_mean[1];

  //real<lower = lambda_mean_prior-lambda_sd_prior,upper=lambda_mean_prior+lambda_sd_prior> lambda_mean[1];
  //vector[numb_year-1] lambda_var; // loop parameter
  real<lower = g_sd_min,upper=1> g_init[1]; //g_mean+g_sd
}

transformed parameters{

  // Declare objects neceseed_sary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N-numb_plot] interaction_effects;
  //vector<lower = 0>[numb_year] g_ei;
    matrix[N-numb_plot,S] alpha_value;
    matrix[N-numb_plot,S] alpha_ei;
  vector[N-numb_plot] seedtotal; // estimated number of seed per year and plot
  vector<lower = 0>[N-numb_plot] lambda_ei; // loop parameter
  matrix[N-numb_plot,S] SpAbundance_scaled;
  //real g_init;
 // implement the biological model
  // for (y in 1:numb_year){
    //g_ei[y] = g_init[1] + beta2[1]*PDSI_mean[y];
  // }
    for(n in 1:N-numb_plot){
    lambda_ei[n] = lambda_mean_prior + beta[1]*PDSI_mean[year[n]];
    //lambda_ei[n] = lambda_mean[1] + lambda_var[year[n]];
      //g_init = g_mean;
  for(s in 1:S){
    SpAbundance_scaled[n,s] = SpAbundance[n,s]*0.6;// 15/25
    alpha_value[n,s] = alpha_initial_mean[s] + (c_mean[s] * (1 - exp(alpha_slope_mean[s] * (SpAbundance_scaled[n,s]-N_opt_mean[s]))))/(1+exp(alpha_slope_mean[s] * (SpAbundance_scaled[n,s]- N_opt_mean[s])));
    alpha_ei[n,s]=alpha_value[n,s]*SpAbundance_scaled[n,s];
    }
   //print(SpAbundance_scaled[n,]);
    //print(alpha_ei[n,]);
    
   interaction_effects[n] = sum(alpha_ei[n,]);
  if (obs_spI[n] > 0){
  // the total number of seed at t+1 is influenced by the number of seed produced at t, which depends up the abundance of neighbours at t, lambda at t
   seedtotal[n] = (1- g_init[1])*seed_s[1]*obs_spI[n]./g_init[1] + obs_spI[n]*lambda_ei[n]*exp(interaction_effects[n]);
    //print(seedtotal[n]);
   }
  else{ // if obsspI at t is 0, yet seedtota can't be equal to 0 if the year after we have a  obsspI at t+1 not equal to 0,
  // all the seeds are in the seed bank and no plants either germinated or survived. . 
  // So we can replace with the median abundance form the previous year
  // the total number of seed at t+1 is influenced by the number of seed produced at t, which depends up the abundance of neighbours at t, lambda at t
  // if the Abundance 
   seedtotal[n] = (1-g_init[1] )*seed_s[1]*SpAbundance_med[year[n],1]./g_init[1]+  SpAbundance_med[year[n],1]*lambda_ei[n]*exp(interaction_effects[n]);
   }
  }

}

model{
  matrix[1,N-numb_plot] stems_spI ; // estimated number of stem per year and plot
  //vector[N-numb_plot] stems_spI;
  seed_s ~  beta(5,5); //normal(seed_s_mean,seed_s_sd); //beta(5,5);
  //beta  ~ normal(beta_prior,1);
  beta  ~ normal(0,1);
  //lambda_mean ~ normal(lambda_mean_prior,lambda_sd_prior);
  g_init ~ beta(5,5); //normal(g_mean,g_sd); // normal distribution with mena 0.5 and variance 0.022  
  //for(y in 1:numb_year-1){
  //lambda_var[year[y]] ~ normal(0,1);
  //}
  // disp_dev ~ normal(0,1);
   for(n in 1:N-numb_plot){
  //PDSI[year[n]] ~ normal(PDSI_mean[year[n]],PDSI_sd[year[n]]); // draw a value of PDSI proper to that year
  stems_spI[1,n] = seedtotal[n]* g_init[1]; // g_ei[year[n]+1] get the number of stem for the year y and plot p for species I 
  // the probability distribtuion of species I for year t follows a poisson distribution include all observation across plots of the year t-1
  //print(stems_spI[n]);
  //if(stems_spI[1,n]> 0){
  obs_spI[n+numb_plot] ~  poisson(stems_spI[1,n]); // draw a value of observed species I for that year 
  //}
  //else{
  // SpAbundance_med[year[n]+1,1] ~ lognormal(stems_spI[n],1/disp_dev[1]^2);
   // }
  }
}
