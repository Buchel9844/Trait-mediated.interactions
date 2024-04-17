// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;  // Number of plots/obs
  int<lower = 1> S;  // Number of species
  int<lower = 1> U;  // Upper bound lambda
  int<lower = 0> Nmax[S]; // density at max fecundity, mean for opt
  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
   int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
}

parameters{
  vector<lower=0,upper =1>[1] lambdas; // intrinsic growth rate - always positive
  vector<lower=0,upper =1>[S] Rho; // density optimal multiplicator
  
  vector<lower=-0.5,upper =0.5>[S] c; //stretching parameters

  vector<lower=-1,upper =1>[S] alpha_initial; // initial effect of j on i - when Nj is minimal
  vector<lower=-0.5,upper =0.5>[S] alpha_slope; // decay - impact of the addition of one individual of j, on the fecundity of i. 

    
  real<lower=0> disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{
 // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  
 // loop parameters
  matrix[N,S] alpha_function_eij;
  matrix[N,S] alpha_value;
  vector[N] lambda_ei;
  vector[S] alpha_slope_ei;
  vector[S] c_ei;
  vector[S] N_opt_i; // Density optimal - always positive

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = U*lambdas[1];
    for(s in 1:S){
    //scaling factor
    alpha_slope_ei[s] = alpha_slope[s] - 0.5; //scaling to have higher values
    c_ei[s] = c[s] - 0.5;
    N_opt_i[s] = Rho[s]*Nmax[s];
    
      alpha_value[i,s]= alpha_initial[s] + (c_ei[s]*(1 - exp( alpha_slope_ei[s]*(SpMatrix[i,s]-N_opt_i[s]))))/(1+exp( alpha_slope_ei[s]*(SpMatrix[i,s]-N_opt_i[s])));
      alpha_function_eij[i,s]= alpha_value[i,s]*SpMatrix[i,s];
      
    }
    
 F_hat[i] = exp(lambda_ei[i] + sum(alpha_function_eij[i,]));

  }

}

model{
  // set regular priors
  alpha_initial ~ normal(0,1);
  alpha_slope ~ normal(0,1);
  lambdas ~ normal(0,1);
  c ~ normal(0, 1);
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  Rho ~ normal(0,1);

 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev^2)^(-1)); 
   }

}

generated quantities{
  vector[N] F_sim;
    if(run_estimation==1){
 for(i in 1:N){
    if(F_hat[i] <= 0) break ;
    F_sim[i] = neg_binomial_2_lpmf(Fecundity[i]|F_hat[i],(disp_dev^2)^(-1));
              }
    }
}
