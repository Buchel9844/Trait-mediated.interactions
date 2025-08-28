Read me file for the article entitled 
"Trait-mediated interactions drive local diversity"
by
Lisa Buche, Oscar Godoy, Lauren Shoemaker, Lauren Hallett, Courtney Taylor, Wing Man Siu, Manuel Sevenello, Peter Vesk, and Margie Mayfield 
Published in 
""

Buche L., Godoy, O. et al (2025), Buchel9844/Facilitation_gradient: Initial release (V.0). Zenodo. DOI:  (version V.1 will be released upon acceptance).
Contact details: Lisa Buche (buchel9844@gmail.com)

Statement of authorship: The research group led by Margaret Mayfield, including Lisa Buche, Manuel Sevenello, Wing Man Siu, and Courtney Taylor, collected the data from Australia. Oscar Godoy and Lisa Buche collected the data from Spain. Lisa Buche built the conceptual idea with input from all co-authors, especially Peter Vesk and Margaret Mayfield. Lisa Buche analyzed the data with substantial input from Peter Vesk, Lauren Shoemaker, and Oscar Godoy. Lisa Buche wrote the first draft with substantial input from Margaret Mayfield, Oscar Godoy, Lauren Shoemaker, Lauren Hallett and Peter Vesk. All co-authors reviewed the manuscript.

Abstract: 


Authorship of data: Members of the Mayfield lab collected the data in Australia. We want to thank John Dwyer, Claire Wainwright, Maia Raymundo, Trace Martyn, Victoria Reynolds, Catherine Bowler, Aubrie James, Abigail Pastore, Manuel Sevenello, and Courtney Taylor for the data collected in the Perenjori region. Oscar Godoy and his lab collected the data in Spain.

Authorship of code: R Code was written by Lisa Buche.


Details of R script and their function in the folder CODE:

1.DataPrep.R - Gathers data from all sources (internal files and online resources for traits). The script is divided between Spain related data gathering and Australia. The scripts creates one list object for each country (i.e., "data/clean.data.aus.RData"), containing all the data separated in different data frame.  


2.NaturalDataFit.R - Fits the bayesian model contained in ("Ricker_fit.stan") stan scripts for each focal species and country. Creates a list object for each fit (i.e., results/Parameters_",Code.focal,"_",country,".Rdata") countianing the fitted paramters and model behavior check. 
      2.1.CheckModelBehavior.R - Extracts the models' convergence, Posterior    predictive distribution of the model fit and explanatory power of the fitted models, which were evaluated according to Root mean squared deviance and the leave-one-out approximation (see Appendix section "Model behavior of species interactions" and Tab S5). 

3.Results.R - Makes some supplementary figures of abundance over time for the appendix (e.g., Fig S1-S3). Gather the parameter's fit into one dataset for each country (i.e.,"results/Param.sigmoid.",country,".csv.gz"). Contains the code to produce some supplementary figures of the estimated intrinsic growth rates of each species over time (Fig S9-S10) and the raw sigmoid (Fig S11-S12).

3.RealisedInteractionResults.R - Compute the interaciton outcome for the low and high neighbor densities (i.e., "results/Theoretical.Int.list.RData"). Creates the network and Fig 2. Also make the summary table of interactions (Tab. S6).

3.2 TraitResults.R - Runs the GLM with BRMS packages to estimated the effect size of the trait variables. It also produce fig 3-4-5 and supplementary figures Fig S13-S14.


