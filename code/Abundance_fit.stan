// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/aseed_sets/case_studies/bayeseed_sparse_regreseed_sion.html

data{
  int<lower = 1> N; // Number of observation total
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int focal_pos; //position of intra in SpAbundance dataset 
  int<lower = 1> numb_year; // number of years
  int<lower = 1> numb_plot; // number of plot observation within each year
 
  matrix[N,S] SpAbundance;  // Abundance of the species in each observation 
   matrix[
  numb_year,1] SpAbundance_med;
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
  real<lower=0,upper=1> seed_s[1];
  real PDSI[numb_year-1];

  real<lower=-1,upper=1> alpha_initial[S];
  real<lower=-1,upper=0> alpha_slope[S];
  real<lower=-1,upper=0> c[S];
  real<lower = 0> N_opt[S];
  
  vector<lower=0>[1] disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  real beta[1];
  real beta2[1];
  real<lower = 0> lambda_init[1];
  real<lower = 0,upper=1> g_init[1];
}

transformed parameters{

  // Declare objects neceseed_sary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N-numb_plot] interaction_effects;
  //  vector<lower = 0>[numb_year] g_ei;
  vector<lower = 0>[numb_year-1] lambda_ei; // loop parameter
  vector<lower = 0>[N-numb_plot] seedtotal; // estimated number of seed per year and plot
  vector[S] alpha_value;
 // implement the biological model
   for (y in 1:numb_year){
    g_ei[y] = g_init[1] + beta2[1]*PDSI_mean[y];
   }
    for(n in 1:N-numb_plot){
    lambda_ei[year[n]] = lambda_init[1] + beta[1]*PDSI_mean[year[n]];
  for(s in 1:S){
    alpha_value[s] = alpha_initial[s] + (c[s] * (1 - exp(alpha_slope[s] * (SpAbundance[n,s]-N_opt[s]))))/(1+exp(alpha_slope[s] * (SpAbundance[n,s] - N_opt[s])));
      }
   interaction_effects[n] = sum(alpha_value * SpAbundance[n,]);
   
  if (obs_spI[n] > 0){
  // the total number of seed at t+1 is influenced by the number of seed produced at t, which depends up the abundance of neighbours at t, lambda at t
   seedtotal[n] = (1- g_ei[year[n])*seed_s[1]*obs_spI[n]/ g_ei[year[n] + obs_spI[n]*lambda_ei[year[n]]*exp(interaction_effects[n]);
    //print(seedtotal[n]);
   }
  else{ // if obsspI at t is 0, yet seedtota can't be equal to 0 if the year after we have a  obsspI at t+1 not equal to 0,
  // all the seeds are in the seed bank and no plants either germinated or survived. . 
  // So we can replace with the median abundance form the previous year
  // the total number of seed at t+1 is influenced by the number of seed produced at t, which depends up the abundance of neighbours at t, lambda at t
  // if the Abundance 
   seedtotal[n] = (1-  g_ei[year[n]])*seed_s[1]*SpAbundance_med[year[n],1]/g_ei[year[n] +  SpAbundance_med[year[n],1]*lambda_ei[year[n]]*exp(interaction_effects[n]);
   }
  }
}

model{
  vector[N-numb_plot] stems_spI; // estimated number of stem per year and plot

   for(s in 1:S){
   N_opt[s] ~ normal(N_opt_mean[1,s], N_opt_sd[1,s]);
   alpha_initial[s] ~ normal(alpha_initial_mean[1,s], alpha_initial_sd[1,s]);  // normal prior on lambda
   alpha_slope[s] ~ normal(alpha_slope_mean[1,s], alpha_slope_sd[1,s]);     // normal prior on alpha
   c[s] ~ normal(c_mean[1,s], c_sd[1,s]);     // normal prior on alpha
     }
    
  disp_dev ~ normal(0,1);
  seed_s ~ beta(5,5); 
  beta  ~ normal(0,1);
  beta2  ~ normal(0,1);
  lambda_init ~ normal(0,1);
  g_init ~ beta(5,5); // normal distribution with mena 0.5 and variance 0.022  

   for(n in 1:N-numb_plot){
  //PDSI[year[n]] ~ normal(PDSI_mean[year[n]],PDSI_sd[year[n]]); // draw a value of PDSI proper to that year
  stems_spI[n] = seedtotal[n]*g_init[1]; // g_ei[year[n]+1] get the number of stem for the year y and plot p for species I 
  // the probability distribtuion of species I for year t follows a poisson distribution include all observation across plots of the year t-1
  //print(stems_spI[n]);
  if(stems_spI[n] > 0)
  obs_spI[n+numb_plot] ~ poisson(stems_spI[n]); // draw a value of observed species I for that year 
  }
}


// generated quantities{
 //  vector[N] stems_spI; // estimated number of stem per year and plot
 //  vector[N] pred_spI;
 //  for(n in 1:N){
  //   stems_spI[n] = seedtotal[n]*g[1];
  // if(stems_spI[n] <= 0) break ;
  //   pred_spI[n] = poisson_lpmf(obs_spI[n]|stems_spI[n]); // draw a value of observed species I for that year 
 //  }
// }