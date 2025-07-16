#Yates Coley
#rycoley@gmail.com
#Updated 2017-7-11
#Annotations updated 12/21/17
#This script will run all scripts for model estimation and preparing patient-level predictions


### WORKFLOW
## 1. Clear workspace
## 2. Define directories, file names
## 3. Load libraries
## 4. Set date of data pull
## 5. Run R scripts for model estimation

### 1. Clear workspace
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
print(args)
#mri_role <- "both" #moderator, both, outcome, 0
mri_role <- args[3]
to_mask <- as.numeric(args[2])
K <- as.numeric(args[1])
workdir <- args[4]
MRI_effect <- args[5] #strength or binary
### 2. Define directories, file names
#workdir <- "/users/zwang3/PAS_11722"
base.location <- workdir
location.of.data <- paste0(base.location, "/data")
location.of.r.scripts <- paste0(base.location, "/code")
location.of.generated.files <- paste0(base.location, "/generated-files")
location.of.generated.folder = paste(location.of.generated.files, "/K=",K,"_MRIeffect=",MRI_effect,sep="")
ifelse(!dir.exists(location.of.generated.folder), dir.create(location.of.generated.folder), FALSE)

# base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/"
# location.of.data <- paste0(base.location, "data")
# location.of.r.scripts <- paste0(base.location, "cleaned_code")
# location.of.generated.files <- paste0(base.location, "generated-files-sh")


name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"
name.of.mri.data <- "MRI_6.15.csv" 

#name.of.pt.data <- "Demographics_10.29.22.csv" 
#name.of.bx.data <- "Biopsy_processed.csv" 
#name.of.psa.data <- "PSA_10.29.22.csv" 
#name.of.tx.data <- "Treatments_10.29.22.csv" 
#name.of.mri.data <- "MRI_10.29.22.csv" 
#name.of.mri.findings.data <- "MRI_findings_10.29.22.csv"
#name.of.targeted.pathology.data <- "Targeted_pathology_10.29.22.csv"



### 3. Load libraries
#### may need to go back and add command for automatic installation - YESS done Zitong
# Package names
packages <- c("lme4",  "dplyr", "tidyr", "readr",
              "splines", "bayesm", "rjags", "R2jags")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
### 5. Run R scripts!

#Load data; tidy, check, and shape
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}

#Load data; tidy, check, and shape
#source("code/data-load-check-and-shaping.R") #need to rerun this with new data
load(paste(location.of.generated.files,"IOP-data-shaping-work-space-6.15-withnewMRI.RData",sep="/"))
to.mask<- to_mask

options(warn = 0)
#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"cv-data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"cv-argument-prep-for-jags.R",sep="/"))

#run model and save results
if(K==1){
  seed <- 2022
}else{
  seed <- to.mask
}

#Define JAGS model
source(paste(location.of.r.scripts,"cv-JAGS-prediction-model.R",sep="/"))

set.seed(seed)
outj<-jags(jags_data, inits=inits, 
           parameters.to.save=params,
           model.file=paste(location.of.r.scripts, 
                            jags_file_name, sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput

if (K > 1){
  for(j in 1:length(out$sims.list)){
    if (names(out$sims.list)[j] == "eta.track") {
      write.csv(out$sims.list[[j]],
                paste(location.of.generated.folder, "/jags-prediction-", names(out$sims.list)[j],"-",mri_role, seed,".csv",sep=""))
    }
  }
}else{
  for(j in 1:length(out$sims.list)){
    #if (names(out$sims.list)[j] == "eta.track") {
      write.csv(out$sims.list[[j]],
                paste(location.of.generated.folder, "/jags-prediction-", names(out$sims.list)[j],"-",mri_role,".csv",sep=""))
    #}
     save(out, file = paste(location.of.generated.folder, "/jags-prediction-alloutputs-", mri_role,".RData",sep=""))
  }
}
