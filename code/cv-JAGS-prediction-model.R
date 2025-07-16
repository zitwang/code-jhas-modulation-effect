if(K>1){
  jags_file_name <- paste0("cv-JAGS-prediction-model-",mri_role,seed,MRI_effect,".txt")
}else{
  jags_file_name <- paste0("JAGS-prediction-model-",mri_role,seed,MRI_effect,".txt")
}


cat(
  "model {

###PRIORS FOR LATENT CLASS MODEL

# priors for sequential model for ordinal eta as outcome
for (index in 1:npred_cancer) {
  cancer_coef_mean[index] ~ dnorm(0,1)
}
cancer_int1 ~ dnorm(-0.16, 1) 
cancer_int2 ~ dnorm(0.71, 1)
cancer_int3 ~ dnorm(1.15, 1)
for(index in 1:npred_cancer) {
  cancer_slope1[index] ~ dnorm(cancer_coef_mean[index], 1)
  cancer_slope2[index] ~ dnorm(cancer_coef_mean[index], 1)
  cancer_slope3[index] ~ dnorm(cancer_coef_mean[index], 1)
}


###PRIORS FOR PSA MIXED MODEL

#model correlated random effects distribution
scale_prior_range <- 100
gaussian_prior_prec <- 0.01
for (index in 1:npred_ranef_psa) {
	for(k in 1:nlevel_cancer_bin) {
		mu_raw[index, k] ~ dnorm(0, gaussian_prior_prec)
	}  
}
	
#same covariance matrix (Sigma_B) across latent classes
Tau_B_raw ~ dwish(I_npred_ranef_psa[,], (npred_ranef_psa+1))
Sigma_B_raw[1:npred_ranef_psa, 1:npred_ranef_psa] <- inverse(Tau_B_raw[1:npred_ranef_psa, 1:npred_ranef_psa])

#residual variance, independent of correlated random effects, same across classes
resid_var_prior_range <- 1
resid_var_psa ~ dunif(0, resid_var_prior_range)
tau_res <- pow(resid_var_psa, -2)

#fixed effects
for(index in 1:npred_fixef_psa) {
	fixef_coefficient[index] ~ dnorm(0, gaussian_prior_prec)
}



##proportional odds regression for biopsy grade
for (index in 1:npred_pgg) {
  pgg_coef_mean[index] ~ dnorm(0,1)
}
for (index in (npred_pgg + 1):(npred_pgg + (nlevel_cancer - 1))) {
  pgg_coef_mean[index] ~ dnorm(0,1)
}
",
if(mri_role %in% c("moderator", "both")){"pgg_coef_mean[npred_pgg + nlevel_cancer] ~ dnorm(0,1)"},

"\n
pgg_int1 ~ dnorm(-0.16, 1) 
pgg_int2 ~ dnorm(0.71, 1)
pgg_int3 ~ dnorm(1.15, 1)
for(index in 1:npred_pgg) {
  pgg_slope1[index] ~ dnorm(pgg_coef_mean[index], 1)
  pgg_slope2[index] ~ dnorm(pgg_coef_mean[index], 1)
  pgg_slope3[index] ~ dnorm(pgg_coef_mean[index], 1)
}
for (index in (npred_pgg + 1):(npred_pgg + (nlevel_cancer - 1))) {
  pgg_slope1[index] ~ dnorm(pgg_coef_mean[index], 1)
  pgg_slope2[index] ~ dnorm(pgg_coef_mean[index], 1)
  pgg_slope3[index] ~ dnorm(pgg_coef_mean[index], 1)
}",

if(mri_role %in% c("moderator", "both")){
  "
  pgg_slope1[npred_pgg + nlevel_cancer] ~ dnorm(pgg_coef_mean[npred_pgg + nlevel_cancer],1)
  pgg_slope2[npred_pgg + nlevel_cancer] ~ dnorm(pgg_coef_mean[npred_pgg + nlevel_cancer],1)
  pgg_slope3[npred_pgg + nlevel_cancer] ~ dnorm(pgg_coef_mean[npred_pgg + nlevel_cancer],1)"
},

"

###LIKELIHOOD
##latent variable for true cancer state (sequential logistic regression)
for(j in 1:npat){ 
  exponent1[j] <- cancer_int1 + inprod(cancer_slope1[1:npred_cancer], modmat_cancer[j,1:npred_cancer])
  exponent2[j] <- cancer_int2 + inprod(cancer_slope2[1:npred_cancer], modmat_cancer[j,1:npred_cancer])
  exponent3[j] <- cancer_int3 + inprod(cancer_slope3[1:npred_cancer], modmat_cancer[j,1:npred_cancer])
  p_eta[j, 1] <- exp(exponent1[j])/(1+exp(exponent1[j]))
  p_eta[j, 2] <- 1/(1+exp(exponent1[j])) * exp(exponent2[j])/(1+exp(exponent2[j]))
  p_eta[j, 3] <- 1/(1+exp(exponent1[j])) * 1/(1+exp(exponent2[j])) * exp(exponent3[j])/(1+exp(exponent3[j]))
  p_eta[j, 4] <- 1- p_eta[j, 1] - p_eta[j, 2] - p_eta[j, 3]
}


for(i in 1:npat_cancer_known) {
  cancer_data[i] ~ dcat(p_eta[i, 1:nlevel_cancer])
  eta[i] <- cancer_data[i]
  eta.bin[i] <- step(eta[i]-2) + 1
  eta_eq2[i] <- ifelse(eta[i] == 2, 1, 0)           ## indicator for eta == 2, ==3, ==4
  eta_eq3[i] <- ifelse(eta[i] == 3, 1, 0)
  eta_eq4[i] <- ifelse(eta[i] == 4, 1, 0)
} #this is for those with path reports from SURG, eta known

for(i in (npat_cancer_known+1):npat) {
  cancer_state[(i-npat_cancer_known)] ~ dcat(p_eta[i, 1:nlevel_cancer])
  eta[i] <- cancer_state[(i-npat_cancer_known)]
  eta.bin[i] <- step(eta[i]-2) + 1
  eta_eq2[i] <- ifelse(eta[i] == 2, 1, 0)           ## indicator for eta == 2, ==3, ==4
  eta_eq3[i] <- ifelse(eta[i] == 3, 1, 0)
  eta_eq4[i] <- ifelse(eta[i] == 4, 1, 0)
}  #for those without SURG

",

if(K>1){
  "eta.track[1:n_mask] <- cancer_state[1:n_mask]"
},

"

##linear mixed effects model for PSA
#generate random intercept and slope for individual given latent class
for (i in 1:npat) {
	B_raw[i, 1:npred_ranef_psa] ~ dmnorm(mu_raw[1:npred_ranef_psa, eta.bin[i]], Tau_B_raw[1:npred_ranef_psa, 1:npred_ranef_psa])
	for(index in 1:npred_ranef_psa) {
	  ranef[i, index] <- B_raw[i, index]
	} 
}

#fit LME
for(j in 1:nobs_psa){
	mu_obs_psa[j] <- inprod(ranef[psa_patient_index_map[j], 1:npred_ranef_psa], modmat_ranef_psa[j, 1:npred_ranef_psa]) + inprod(fixef_coefficient[1:npred_fixef_psa], modmat_fixef_psa[j,1:npred_fixef_psa])
	log_psa_data[j] ~ dnorm(mu_obs_psa[j], tau_res) 
}

### BIOPSY OUTCOMES AND SURGERY RECEIVED

##proportional odds regression for biopsy upgrading
for(j in 1:npat_pgg){ 
  pgg_exp1[j] <- pgg_int1 + inprod(pgg_slope1[1:npred_pgg], modmat_pgg[j,1:npred_pgg])+
	               pgg_slope1[(npred_pgg + 1)] * step(eta[pgg_patient_index_map[j]] - 2) + 
	               pgg_slope1[(npred_pgg + 2)] * step(eta[pgg_patient_index_map[j]] - 3) + 
	               pgg_slope1[(npred_pgg + 3)] * step(eta[pgg_patient_index_map[j]] - 4)",
if(mri_role %in% c("moderator", "both")){
  if(MRI_effect == "binary"){
    "+ pgg_slope1[(npred_pgg + 4)] * (step(eta[pgg_patient_index_map[j]] - 2)-(1-step(eta[pgg_patient_index_map[j]] - 2))) * step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
    ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  }else{
    "+ pgg_slope1[(npred_pgg + 4)] * (pgg_pirads_data_m2[pgg_patient_index_map[j]]) * (step(eta[pgg_patient_index_map[j]] - 2)-(1-step(eta[pgg_patient_index_map[j]] - 2))) * step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
    ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  }
},
"
  pgg_exp2[j] <- pgg_int2 + inprod(pgg_slope2[1:npred_pgg], modmat_pgg[j,1:npred_pgg])+
	               pgg_slope2[(npred_pgg + 1)] * step(eta[pgg_patient_index_map[j]] - 2) + 
	               pgg_slope2[(npred_pgg + 2)] * step(eta[pgg_patient_index_map[j]] - 3) + 
	               pgg_slope2[(npred_pgg + 3)] * step(eta[pgg_patient_index_map[j]] - 4)",
if(mri_role %in% c("moderator", "both")){
  if(MRI_effect == "binary"){
    "+ pgg_slope2[(npred_pgg + 4)] * (pgg_pirads_data_m2[pgg_patient_index_map[j]]) * (step(eta[pgg_patient_index_map[j]] - 3)-(1-step(eta[pgg_patient_index_map[j]] - 3))) * step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
	     ## logitP(pgg = 2|pgg >=2) w/ (pirads - 2)*[1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
  }else{
    "+ pgg_slope2[(npred_pgg + 4)] *  (step(eta[pgg_patient_index_map[j]] - 3)-(1-step(eta[pgg_patient_index_map[j]] - 3))) * step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
	     ## logitP(pgg = 2|pgg >=2) w/ [1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
  }
},
"
 pgg_exp3[j] <- pgg_int3 + inprod(pgg_slope3[1:npred_pgg], modmat_pgg[j,1:npred_pgg])+
	               pgg_slope3[(npred_pgg + 1)] * step(eta[pgg_patient_index_map[j]] - 2) + 
	               pgg_slope3[(npred_pgg + 2)] * step(eta[pgg_patient_index_map[j]] - 3) + 
	               pgg_slope3[(npred_pgg + 3)] * step(eta[pgg_patient_index_map[j]] - 4)",
if(mri_role %in% c("moderator", "both")){
  if(MRI_effect == "binary"){
    "+ pgg_slope3[(npred_pgg + 4)] *  (step(eta[pgg_patient_index_map[j]] - 4)-(1-step(eta[pgg_patient_index_map[j]] - 4))) * 
    step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
	   ## logitP(pgg =3|pgg>=3) w/ [1(eta >= 4)-1(eta<=3)]*1(pirads >= 3)"
  }else{
    "+ pgg_slope3[(npred_pgg + 4)] * (pgg_pirads_data_m2[pgg_patient_index_map[j]]) * (step(eta[pgg_patient_index_map[j]] - 4)-(1-step(eta[pgg_patient_index_map[j]] - 4))) * 
    step(pgg_pirads_data[pgg_patient_index_map[j]] - 3)  
	   ## logitP(pgg =3|pgg>=3) w/ (pirads - 2)*[1(eta >= 4)-1(eta<=3)]*1(pirads >= 3)"
  }
},
"
  
  p_pgg[j, 1] <- exp(pgg_exp1[j])/(1+exp(pgg_exp1[j]))
  p_pgg[j, 2] <- 1/(1+exp(pgg_exp1[j])) * exp(pgg_exp2[j])/(1+exp(pgg_exp2[j]))
  p_pgg[j, 3] <- 1/(1+exp(pgg_exp1[j])) * 1/(1+exp(pgg_exp2[j])) * exp(pgg_exp3[j])/(1+exp(pgg_exp3[j]))
  p_pgg[j, 4] <- 1- p_pgg[j, 1] - p_pgg[j, 2] - p_pgg[j, 3]
}

for(i in 1:npat_pgg) {
  pgg_data[i] ~ dcat(p_pgg[i,1:nlevel_cancer])
}",

if(mri_role %in% c("outcome", "both")){
  "#Prior for mri pirads data
for (index in 1:npred_pirads) {
  pirads_coef_mean[index] ~ dnorm(0,1)
}
for (index in (npred_pirads + 1):(npred_pirads + (nlevel_cancer - 1))) {
  pirads_coef_mean[index] ~ dnorm(0,1)
}
pirads_int1 ~ dnorm(-0.69, 1) 
pirads_int2 ~ dnorm(-0.51, 1)

for(index in (npred_pirads + 1):(npred_pirads + (nlevel_cancer - 1))) {
  pirads_slope1[index] ~ dnorm(pirads_coef_mean[index], 1)
  pirads_slope2[index] ~ dnorm(pirads_coef_mean[index], 1)
}

### SEQUENTIAL LOGISTIC MODEL FOR MRI PI-RADS SCORE
for(j in 1:npat_pirads){ 
  pirads_exp1[j] <- pirads_int1 +
	               pirads_slope1[(npred_pirads + 1)] * step(eta[pirads_patient_index_map[j]] - 2) + 
	               pirads_slope1[(npred_pirads + 2)] * step(eta[pirads_patient_index_map[j]] - 3) + 
	               pirads_slope1[(npred_pirads + 3)] * step(eta[pirads_patient_index_map[j]] - 4)
  pirads_exp2[j] <- pirads_int2 +
	               pirads_slope2[(npred_pirads + 1)] * step(eta[pirads_patient_index_map[j]] - 2) + 
	               pirads_slope2[(npred_pirads + 2)] * step(eta[pirads_patient_index_map[j]] - 3) + 
	               pirads_slope2[(npred_pirads + 3)] * step(eta[pirads_patient_index_map[j]] - 4)
  p_pirads[j, 1] <- exp(pirads_exp1[j])/(1+exp(pirads_exp1[j]))
  p_pirads[j, 2] <- 1/(1+exp(pirads_exp1[j])) * exp(pirads_exp2[j])/(1+exp(pirads_exp2[j]))
  p_pirads[j, 3] <- 1- p_pirads[j, 1] - p_pirads[j, 2] 
}

for(i in 1:npat_pirads) {
  pirads_data[i] ~ dcat(p_pirads[i,1:(nlevel_cancer-1)])
  }"
},

"}",
  file = paste(location.of.r.scripts, 
                jags_file_name, sep="/"))


