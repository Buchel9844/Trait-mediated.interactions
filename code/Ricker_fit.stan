// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of observations within each year
  int<lower = 1> S; // Number of plant species (same across years for consistency)
  int<lower = 1> year[N]; // Indicator variable for the year each observations
  int<lower = 1> Y; // number of years
  int<lower = 1> lambda_max; // max fecundity expected
  int<lower = 0> lambda_min; // min fecundity expected
  real N_opt_prior[S];
  int Fecundity[N];  // Fecundity of the focal species in each observation
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (nncluding abundances of non-focal individuals of the focal species)
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations

}

parameters{
  vector<lower=lambda_min,upper=lambda_max>[1] lambda_mean;
  vector[Y] lambda_sd;

  vector<lower=-1,upper =0>[S] c; //stretching parameters

  vector<lower=-1,upper =0>[S] alpha_slope; // decay - impact of the addition of one individual of j, on the fecundity of i. 

  vector<lower=-1,upper =1>[S] alpha_initial; // initial effect of j on i - when Nj is minimal
  vector<upper=0>[1] alpha_init_intra;
  vector<lower=0>[S] N_opt_i; 
  
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
  vector<lower=lambda_min,upper=lambda_max>[N] lambda_ei;
  matrix[N,S] alpha_value;
  

 // implement the biological model

  for(n in 1:N){
    lambda_ei[n] = (lambda_mean[1] + lambda_sd[year[n]]);

    for(s in 1:S){
    alpha_value[n,s] = (1-Intra[s])*alpha_initial[s] + (Intra[s])*alpha_init_intra[1] + (c[s]*(1 - exp( alpha_slope[s]*(SpMatrix[n,s]-N_opt_i[s]))))/(1+exp(alpha_slope[s]*(SpMatrix[n,s]-N_opt_i[s])));
    }
    
    interaction_effects[n] = sum(alpha_value[n,] .* SpMatrix[n,]);
     
    F_hat[n] =  lambda_ei[n]*exp(interaction_effects[n]);
   }
}

model{
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  for(s in 1:S){
  N_opt_i[s]~ normal(N_opt_prior[s],1); // for definition of neutral density
  }
  
  alpha_initial ~ normal(0,0.1);
  alpha_init_intra ~ normal(-0.1,0.1);
  alpha_slope ~ normal(-0.2,0.2);
  c ~ normal(0,0.1);

  lambda_mean ~ cauchy(0,10); //normal(lambda_prior, lambda_prior_sd);
    
   for(y in 1:Y){
     lambda_sd[y] ~  normal(0,10); //normal(0, lambda_prior_sd);
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
