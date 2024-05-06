// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of observations within each year
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int<lower = 1> year[N]; // Indicator variable for the year each observations
  int<lower = 1> Y; // number of years

  int Fecundity[N];  // Fecundity of the focal species in each observation
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (nncluding abundances of non-focal individuals of the focal species)
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations

  int Inclusion_alpha_initial[Y,S];  // Boolean indicator variables to identify the plant species
  int Inclusion_alpha_slope[Y,S];  // Boolean indicator variables to identify the plant species
  int Inclusion_c[Y,S];  // Boolean indicator variables to identify the plant species
}

parameters{
  real<lower=0> lambda_mean[1];
  real lambda_sd[Y];

  real<lower=-1,upper=1> alpha_initial[S];
  real<lower=-1,upper=1> alpha_initial_hat[Y,S]; 
  
  real<lower=-1,upper=0> alpha_slope[S];
  real<lower=-1,upper=0> alpha_slope_hat[Y,S]; 
  
  real<lower=-1,upper=0> c[S];
  real<lower=-1,upper=0> c_hat[Y,S]; 
  
  real<lower = 0> N_opt_mean[S];
  
    vector<lower=0>[1] disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{

  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (nnteraction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  //vector[N] pollinator_effects;
  
  // loop parameters
  vector[N] lambda_ei;
  matrix[N,S] alpha_initial_ei;
  matrix[N,S] alpha_slope_ei;
  matrix[N,S] c_ei;
  matrix[N,S] N_opt_i;
  matrix[N,S] alpha_value;
  

 // implement the biological model

  for(n in 1:N){
    lambda_ei[n] = lambda_mean[1] + lambda_sd[year[n]];

    for(s in 1:S){
      
    alpha_initial_ei[n,s] = alpha_initial[s] + Inclusion_alpha_initial[year[n],s]*alpha_initial_hat[year[n],s];
    
    alpha_slope_ei[n,s] =  alpha_slope[s] + Inclusion_alpha_slope[year[n],s]*alpha_slope_hat[year[n],s];
    
    c_ei[n,s] =  c[s] + Inclusion_c[year[n],s]*c_hat[year[n],s];
    
    N_opt_i[n,s] = N_opt_mean[s];
    
    alpha_value[n,s] = alpha_initial_ei[n,s] + (c_ei[n,s]*(1 - exp( alpha_slope_ei[n,s]*(SpMatrix[n,s]-N_opt_i[n,s]))))/(1+exp(alpha_slope_ei[n,s]*(SpMatrix[n,s]-N_opt_i[n,s])));
    }
    
    interaction_effects[n] = sum(alpha_value[n,] .* SpMatrix[n,]);
     
    F_hat[n] =  exp(lambda_ei[n] + interaction_effects[n]);
   }
}

model{
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  N_opt_mean ~ normal(0,1); // for definition of neutral density
  
  alpha_initial ~ normal(0,1);
  alpha_slope ~ normal(0,1);
  c ~ normal(0,1);

    lambda_mean ~ normal(0, 1);
    
   for(y in 1:Y){
     lambda_sd[y] ~ normal(0, 1);
      
    alpha_initial_hat[y,] ~ normal(0,1);

    alpha_slope_hat[y,] ~ normal(0,1);

    c_hat[y,] ~ normal(0,1);
    }
    
 for(n in 1:N){
  Fecundity[n] ~ neg_binomial_2(F_hat[n],(disp_dev[1]^2)^(-1)); 
   }
}

generated quantities{
  vector[N] F_sim;
 for(i in 1:N){
    if(F_hat[i] <= 0) break ;
    F_sim[i] = neg_binomial_2_lpmf(Fecundity[i]|F_hat[i],(disp_dev[1]^2)^(-1));
              }
}
