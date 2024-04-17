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
  int<lower = 1> U;  // Upper bound lambda
  int<lower = 0> Nmax[S];  // vector of "optimal" neighbours density for fecundity of focal's
  
  // The below values define the regularized horseshoe priors used for species-specific parameters
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values

}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  real<lower=0,upper=1> lambdas[Y];

  real<lower=-1,upper=1> alpha_initial_intra[Y]; // intra- initial value = density-indep parameter
    real<lower=-1,upper=1> alpha_slope_intra[Y]; // intra- slope = density-dep parameter
      real<lower=-1,upper=1> c_intra[Y]; // intra- stretching parameter
  
  real<lower=-1,upper=1> alpha_initial[S];
  matrix[Y,S] alpha_initial_hat_tilde; // direct interaction inter plants - plants ; species -specific term 
  matrix<lower = 0>[Y,S] initial_hat_shrinkage; // direct interaction inter plants - plants, shrinkage
  
  real<lower=-1,upper=1> alpha_slope[S];
  matrix[Y,S] alpha_slope_hat_tilde; // direct interaction inter plants - plants ; species -specific term 
  matrix<lower = 0>[Y,S] slope_hat_shrinkage; // direct interaction inter plants - plants, shrinkage
 
  real<lower=-1,upper=1> c[S];
  matrix[Y,S] c_hat_tilde; // direct interaction inter plants - plants ; species -specific term 
  matrix<lower = 0>[Y,S] c_hat_shrinkage; // direct interaction inter plants - plants, shrinkage
 
  real<lower = 0> Rho[S];

  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
  
  real<lower=0> disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  
  matrix[Y,S] alpha_initial_hat; // direct interaction inter plants
  matrix[Y,S] initial_hat_shrinkage_tilde; // direct interaction inter plants
  
  matrix[Y,S] alpha_slope_hat; // direct interaction inter plants
  matrix[Y,S] slope_hat_shrinkage_tilde; // direct interaction inter plants
  
  matrix[Y,S] c_hat; // direct interaction inter plants
  matrix[Y,S] c_hat_shrinkage_tilde; // direct interaction inter plants
  
  
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
  
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // Compute prior of sparsity approach - this calculation follows equation 2.8 in Piironen and Vehtari 2013
    for(y in 1:Y){
    for(s in 1:S){ // direct interaction INTER plants
      initial_hat_shrinkage_tilde[y,s] = sqrt( c2 * square(initial_hat_shrinkage[y,s]) / (c2 + square(tau) * square(initial_hat_shrinkage[y,s])) );
      alpha_initial_hat[y,s] = tau * initial_hat_shrinkage_tilde[y,s] * alpha_initial_hat_tilde[y,s];
      
      slope_hat_shrinkage_tilde[y,s] = sqrt( c2 * square(slope_hat_shrinkage[y,s]) / (c2 + square(tau) * square(slope_hat_shrinkage[y,s])) );
      alpha_slope_hat[y,s] = tau * slope_hat_shrinkage_tilde[y,s] * alpha_slope_hat_tilde[y,s];
      
      c_hat_shrinkage_tilde[y,s] = sqrt( c2 * square(c_hat_shrinkage[y,s]) / (c2 + square(tau) * square(c_hat_shrinkage[y,s])) );
      c_hat[y,s] = tau * c_hat_shrinkage_tilde[y,s] * c_hat_tilde[y,s];
    
    }
    }
      
 // implement the biological model

  for(n in 1:N){
         lambda_ei[n] = U*lambdas[year[n]];
    for(s in 1:S){
    
    alpha_initial_ei[n,s] = (1-Intra[s]) * alpha_initial[s] + Intra[s] * alpha_initial_intra[year[n]] + (1-Intra[s])* alpha_initial_hat[year[n],s];
    
    alpha_slope_ei[n,s] = (1-Intra[s]) * alpha_slope[s] + Intra[s] * alpha_slope_intra[year[n]] + (1-Intra[s])* alpha_slope_hat[year[n],s];
    
    c_ei[n,s] = (1-Intra[s]) * c[s] + Intra[s] * c_intra[year[n]] + (1-Intra[s])* c_hat[year[n],s];
    
    N_opt_i[n,s] = Rho[s]*Nmax[s];
    
    alpha_value[n,s] = alpha_initial_ei[n,s] + (c_ei[n,s]*(1 - exp( alpha_slope_ei[n,s]*(SpMatrix[n,s]-N_opt_i[n,s]))))/(1+exp(alpha_slope_ei[n,s]*(SpMatrix[n,s]-N_opt_i[n,s])));
    }
    
    interaction_effects[n] = sum(alpha_value[n,] .* SpMatrix[n,]);
     
    F_hat[n] =  exp(lambda_ei[n] + interaction_effects[n]);
   }
}

model{
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  Rho ~ normal(0,1); // for definition of neutral density
  
  alpha_initial ~ normal(0,1);
  alpha_slope ~ normal(0,1);
  c ~ normal(0,1);
  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

    tau_tilde ~ cauchy(0,1);
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
    
   for(y in 1:Y){
    alpha_initial_intra[y] ~ normal(0, 1); 
    alpha_slope_intra[y] ~ normal(0, 1);
    c_intra[y] ~ normal(0, 1);
      
    lambdas[y] ~ normal(0, 1);
    alpha_initial_hat_tilde[y,] ~ normal(0,1);
    initial_hat_shrinkage[y,] ~ cauchy(0,1);

    alpha_slope_hat_tilde[y,] ~ normal(0,1);
    slope_hat_shrinkage[y,] ~ cauchy(0,1);
     
     c_hat_tilde[y,] ~ normal(0,1);
     c_hat_shrinkage[y,] ~ cauchy(0,1);
  }
    
 for(n in 1:N){
  Fecundity[n] ~ neg_binomial_2(F_hat[n],(disp_dev^2)^(-1)); 
   }
}
