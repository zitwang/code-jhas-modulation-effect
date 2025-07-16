#### WORKFLOW
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked
### 4. Define other jags settings

library(bayesm) #to define rwishart()

### 1. Define data to be sent to jags function ----------

#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4
nlevel_cancer_bin <- 2

jags_data_init <-list(
  ## cancer model
  nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=npat,
  cancer_data=cancer_data, npat_cancer_known=npat_cancer_known,
  modmat_cancer=modmat_cancer, npred_cancer=npred_cancer,
  n_mask=n_mask,
  ## psa model
  nobs_psa=nobs_psa, log_psa_data=log_psa_data, psa_patient_index_map=psa_patient_index_map,
  modmat_ranef_psa=modmat_ranef_psa, modmat_fixef_psa=modmat_fixef_psa, 
  npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
  I_npred_ranef_psa=diag(npred_ranef_psa),
  ## pgg model
  pgg_data=pgg_data, npat_pgg=npat_pgg, pgg_patient_index_map=pgg_patient_index_map,
  modmat_pgg=modmat_pgg, npred_pgg=npred_pgg) 

if(mri_role == 0){
  jags_data <- jags_data_init
}else{
  if(mri_role == "moderator"){
    jags_data_append <- list(pgg_pirads_data = pgg_pirads_data,
                             pgg_pirads_data_m2 = pgg_pirads_data_m2)
  }else if(mri_role == "outcome"){
    jags_data_append <- list(pirads_data = pirads_data, npat_pirads = npat_pirads,
                             pirads_patient_index_map=pirads_patient_index_map,
                             npred_pirads = npred_pirads)
  }else(
    jags_data_append <- list(pgg_pirads_data = pgg_pirads_data,
                             pgg_pirads_data_m2 = pgg_pirads_data_m2,
                             pirads_data = pirads_data, npat_pirads = npat_pirads,
                             pirads_patient_index_map=pirads_patient_index_map,
                             npred_pirads = npred_pirads)
  )
  jags_data <- append(jags_data_init, jags_data_append)
}





### 2. Initialize model parameters ---------------
inits <- function(){
  cancer_state <- pt.data$bx.pgg[is.na(cancer_data)]
  #latent cancer model
  cancer_int1 <- cancer_int2 <-cancer_int3 <-rnorm(1,0,1)
  cancer_slope1 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope2 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope3 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_coef_mean = rnorm(npred_cancer,mean=0,sd=0.25)
  
  #psa model
  scale_ranef_mean_psa <- rlnorm(npred_ranef_psa)
  mu_raw <- as.matrix(cbind(rnorm(npred_ranef_psa),rnorm(npred_ranef_psa)))
  Tau_B_raw <- rwishart((npred_ranef_psa+1),diag(npred_ranef_psa)*var_vec)$W
  resid_var_psa <- min(rlnorm(1),1)
  fixef_coefficient <- rnorm(npred_fixef_psa)
  
  #bgg results model
  pgg_int1 <- pgg_int2 <-pgg_int3 <-rnorm(1,0,1)
  if(mri_role==0 | mri_role =="outcome"){
    pgg_coef_mean = rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
    pgg_slope1 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
    pgg_slope2 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
    pgg_slope3 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
  }else{
    pgg_coef_mean = rnorm(npred_pgg + (nlevel_cancer),mean=0,sd=0.25)
    pgg_slope1 <- rnorm(npred_pgg + (nlevel_cancer),mean=0,sd=0.25)
    pgg_slope2 <- rnorm(npred_pgg + (nlevel_cancer),mean=0,sd=0.25)
    pgg_slope3 <- rnorm(npred_pgg + (nlevel_cancer),mean=0,sd=0.25)
  }
  
  #mri as outcome model 
  if(mri_role %in% c("outcome", "both")){
    pirads_int1 <- pirads_int2 <- rnorm(1,0,1)
    pirads_slope1 <- rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25)
    pirads_slope2 <- rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25)
    pirads_coef_mean = rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25) 
  }
  
  inits_list <- list(cancer_state=cancer_state,
                     cancer_int1=cancer_int1, cancer_int2=cancer_int2, cancer_int3=cancer_int3, 
                     cancer_slope1=cancer_slope1,cancer_slope2=cancer_slope2,cancer_slope3=cancer_slope3,
                     cancer_coef_mean = cancer_coef_mean,
                     
                     scale_ranef_mean_psa=scale_ranef_mean_psa, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, 
                     resid_var_psa=resid_var_psa, fixef_coefficient=fixef_coefficient,
                     
                     pgg_int1=pgg_int1, pgg_int2=pgg_int2, pgg_int3=pgg_int3, 
                     pgg_slope1=pgg_slope1,pgg_slope2=pgg_slope2,pgg_slope3=pgg_slope3,
                     pgg_coef_mean = pgg_coef_mean)
  
  if(mri_role %in% c("outcome", "both")){
    init_list_append <- list(pirads_int1 = pirads_int1, pirads_int2= pirads_int2,
                             pirads_slope1=pirads_slope1,pirads_slope2=pirads_slope2,
                             pirads_coef_mean = pirads_coef_mean)
  }
  
  if(mri_role %in% c("outcome", "both")){
    append(inits_list, init_list_append)
  }else{
    inits_list
  }
}

### 3. Define parameters to be tracked ---------------
params <- c("eta.track",
            "cancer_int1", "cancer_int2", "cancer_int3",
            "cancer_slope1", "cancer_slope2", "cancer_slope3", 
            "cancer_state",
            
            "scale_ranef_mean_psa", "ranef_intercept", "ranef_slope",
            "ranef_var_intercept", "ranef_var_slope", "ranef_cov","resid_var_psa",  
            "ranef",
            "fixef_coefficient",
            
            "pgg_int1", "pgg_int2", "pgg_int3",
            "pgg_slope1", "pgg_slope2", "pgg_slope3"
)

if(mri_role %in% c("outcome", "both")){
  newparams <- c("pirads_int1", "pirads_int2", "pirads_int3",
                 "pirads_slope1", "pirads_slope2", "pirads_slope3" )
  params <- c(params, newparams)
}



### 4. Define other jags settings -------------------- 
#### it is possible for chain to converge with many fewer iterations. Check the appropriate number of iterations for your dataset.
#### change length; burn-in; number thinned; number of chains
## n.iter 10k minimum
## n.burnin 10k minimum
## n.chains == # cores - 1

# change length; burn-in; number thinned; number of chains
n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1






