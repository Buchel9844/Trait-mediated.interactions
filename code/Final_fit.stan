// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int RemoveH; // Remove Herbivor
  int RemoveFV; // Remove FloralVisitors
  int RemoveFvH; // Remove FloralVisitors

  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  int<lower = 1> U;  // Upper bound lambda
  
  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  matrix[N,H] SpMatrix_H; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_FV; // Matrix of abundances for each floral visitors species
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations

  
  matrix[S,S] matrix_HOIs_plant[N]; // Matrix of abundances for each plant species with each other plant species
  matrix[S,H] matrix_HOIs_ijh[N]; // Matrix of abundances for each herbivores species and competitor
  matrix[S,FV] matrix_HOIs_ijf[N]; // Matrix of abundances for each herbivores species and competitor

   // vector telling which interactions to include
   
  int Inclusion_ij[1,S];  // Boolean indicator variables to identify the plant species
  int Inclusion_FV[1,FV]; // Boolean indicator variables to identify the floral visitor species
  int Inclusion_H[1,H]; // Boolean indicator variables to identify the herbivore species
  int beta_Inclusion_plant[S,S];  // Boolean indicator variables to identify the HOIs plant-plant
  int beta_Inclusion_FV[S,FV];  // Boolean indicator variables to identify the HOIs 2 plants-FV
  int beta_Inclusion_H[S,H];  // Boolean indicator variables to identify the HOIs 2 plants- H

}

parameters{
  vector<lower=0,upper =1>[1] lambdas;
  vector<lower=-1,upper=1>[1] alpha_generic;
  vector<lower=-1,upper=1>[1] alpha_intra;
  vector<lower=-1,upper=1>[S] alpha_hat_ij;
  
  matrix<lower=-1,upper=1>[S,S] beta_plant_hat_ijk; // HOIs plants

  real gamma_H_generic[1]; // direct interaction plants - herbivores ; generic
  vector<lower=-1,upper=1>[H] gamma_H_hat_ih; // direct interaction plants - herbivores ; species -specific term

  real gamma_FV_generic[1]; // direct interaction plants - FV ; generic
  vector<lower=-1,upper=1>[FV] gamma_FV_hat_if; // direct interaction plants - FV ; species -specific term

  matrix<lower=-1,upper=1>[S,H] beta_H_hat_ijh; // HOIs herbivores 

  matrix<lower=-1,upper=1>[S,FV] beta_FV_hat_ijf; // HOIs floral visitors
  
  vector<lower=0>[1] disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)

}

transformed parameters{
  
 // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  vector[N] HOI_effects;
  //vector[N] pollinator_effects;
  
  
  // loop parameters
  vector[N] lambda_ei;
  matrix[N,S] alpha_eij;
  matrix[N,H] gamma_H_eih; //gamma herbivores of focal i
  matrix[N,FV] gamma_FV_eif; //gamma floral visitors of focal i
  

    // the matrix of interaction of plant-plant-plant HOIs was not created prior to the script
  matrix[S,S] beta_ijk; // HOIs plant 
    matrix[N,S] matrix_beta_ijk;
  matrix[S,FV] beta_F_ijf;
    matrix[N,S] matrix_beta_F_ijf;
  matrix[S,H] beta_H_ijh; 
    matrix[N,S] matrix_beta_H_ijh;
  
 // implement the biological model based on the different scenarios
 
  for(i in 1:N){
    //Common paramneters in all scenarios - that is regardless of pollinator and/or herbivores presence or absence
     lambda_ei[i] = U*lambdas[1];
        for(s in 1:S){
      alpha_eij[i,s] = (1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) *alpha_hat_ij[s]*Inclusion_ij[1,s];
      
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = plant species  
        beta_ijk[s,k] = beta_Inclusion_plant[s,k]* beta_plant_hat_ijk[s,k];
        }
        matrix_beta_ijk[i,s] = sum(beta_ijk[s,].* matrix_HOIs_plant[i,s]);
   
        }
       
    if(RemoveFvH == 1){
    HOI_effects[i] = sum(matrix_beta_ijk[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
     
    }else{
    if(RemoveFV ==1){
      for(s in 1:S){
          for(h in 1:H){// for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = plant species
        beta_H_ijh[s,h] =  beta_Inclusion_H[s,h]*beta_H_hat_ijh[s,h];
          }
        matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,].* matrix_HOIs_ijh[i,s]);
      }
     for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H_generic[1] + Inclusion_H[1,h]*gamma_H_hat_ih[h];
      }
    HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_H_ijh[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) + sum(gamma_H_eih[i,] .* SpMatrix_H[i,]);
    }
    if(RemoveH ==1){
      for(s in 1:S){
         for(fv in 1:FV){
        beta_F_ijf[s,fv] = beta_Inclusion_FV[s,fv]*beta_FV_hat_ijf[s,fv];
          }
        matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,].* matrix_HOIs_ijf[i,s]);
      }
      for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV_generic[1] + Inclusion_FV[1,fv]*gamma_FV_hat_if[fv];
      }
        
    HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_F_ijf[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) + sum( gamma_FV_eif[i,] .* SpMatrix_FV[i,]);
       }
    }  
    if(RemoveFvH ==0){
      for(s in 1:S){
        for(h in 1:H){// for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = plant species
        beta_H_ijh[s,h] = beta_Inclusion_H[s,h]*beta_H_hat_ijh[s,h];
        }
        matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,].* matrix_HOIs_ijh[i,s]);
   
        for(fv in 1:FV){
        beta_F_ijf[s,fv] =  beta_Inclusion_FV[s,fv]*beta_FV_hat_ijf[s,fv];
        }
        matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,].* matrix_HOIs_ijf[i,s]);
      }
      for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H_generic[1] + Inclusion_H[1,h]*gamma_H_hat_ih[h];
      }
      for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV_generic[1] + Inclusion_FV[1,fv]*gamma_FV_hat_if[fv];
      }
       
    HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_F_ijf[i,]) +  sum(matrix_beta_H_ijh[i,]) ;
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) +  sum(gamma_H_eih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_eif[i,] .* SpMatrix_FV[i,]);
    }
 
          F_hat[i] = exp(lambda_ei[i] + interaction_effects[i] + HOI_effects[i]);

  }

    
}

model{
  // set regular priors
  alpha_generic ~ normal(0,0.1);
  alpha_intra ~ normal(0,0.1);
  gamma_FV_generic ~ normal(0,0.1);
  gamma_H_generic ~ normal(0,0.1);
    
  lambdas ~ normal(0, 1);
  alpha_hat_ij ~ normal(0,0.1);
  gamma_FV_hat_if ~ normal(0,0.1);
  gamma_H_hat_ih ~ normal(0,0.1);
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi

  for(s in 1:S){
  beta_plant_hat_ijk[,s] ~ normal(0,0.1); // the ensemble of the species-specific beta related to one generic beta follow a normal distribution
 
  beta_FV_hat_ijf[s,] ~ normal(0,0.1);
    
    beta_H_hat_ijh[s,] ~ normal(0,0.1);
  }

    
 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev[1]^2)^(-1)); 
   }
}
