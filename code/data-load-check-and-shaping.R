#Yates Coley
#rycoley@gmail.com
#LAST UPDATED and annotations updated 12/21/17
#This script will load .csv files, check values, and shape data for model estimation



### HASHMARK EXPLANATION
# comment explaining code
## workflow description
### header for section of code
#### comment for user
##### comment for yates to consider later



### WORKFLOW
## 1. Load all data files. Rename unique id
## 2. Start pt.data file. Include risk class, dob, date of death.
## 3. Identify diagnostic biopsy for all patients. Put bx dx date into pt.data. Limit analysis to patients with GS<7 at dx.
## 4. Treatment data shaping
## 5. PSA data shaping
## 6. Biopsy data shaping
## 7. Calculate total follow-up time. Categorize patients based on outcomes, LTF
## 8. Order patients based on observed eta
## 9. Additional Shaping biopsy data. Put data into annual intervals
## 10. Laterality data
## 11. Save datasets.



### Data check function
options(warn=1)
data.check <- function(condition, message){
	if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
	stopifnot(condition)}





### 1. Load all data files. Rename unique id. ----------
demo.data <- read.csv(paste(location.of.data, name.of.pt.data, sep="/"))
bx.data <- read.csv(paste(location.of.data, name.of.bx.data, sep="/"))
psa.data <- read.csv(paste(location.of.data, name.of.psa.data, sep="/"))
tx.data <- read.csv(paste(location.of.data, name.of.tx.data, sep="/"))

#number of patients
n<-dim(demo.data)[1]

#Rename the uniuqe pt id just for brevity
demo.data$id<-demo.data$clinical_PTnum
bx.data$id<-bx.data$clinical_PTnum
psa.data$id<-psa.data$clinical_PTnum
tx.data$id<-tx.data$clinical_PTnum

## Basic checks on files
data.check(condition=as.logical(n==length(unique(demo.data$id))), message="Duplicated PTnum in patient-level data frame.")
##### Other checks needed?

#limit to patients with non-missing clinical_PTnum
demo.data<-demo.data[!is.na(demo.data$clinical_PTnum),]

#number of patients
(n<-dim(demo.data)[1])

data.check(condition=as.logical(0==sum(is.na(demo.data$clinical_PTnum))),
           message="Missing clinical_PTnum in demographic data")
data.check(condition=as.logical(0==sum(is.na(bx.data$clinical_PTnum))),
           message="Missing clinical_PTnum in biopsy data")
data.check(condition=as.logical(0==sum(is.na(psa.data$clinical_PTnum))),
           message="Missing clinical_PTnum in PSA data")
data.check(condition=as.logical(0==sum(is.na(tx.data$clinical_PTnum))),
           message="Missing clinical_PTnum in treatment data")



### 2. Start pt.data file. Include dob, date of death.--------
pt.data <- as.data.frame(demo.data$id)
names(pt.data) <- "id"
pt.data$clinical_PTnum <- demo.data$clinical_PTnum #want to keep this for referencing patients later



#PSAD criteria (PSAD<=0.15)
#data.check(condition=as.logical(0==sum(is.na(demo.data$PSAD_criteria))), message="Missing PSAD criteria values.")
pt.data$lr.psad<-as.numeric(demo.data$PSAD_criteria==0)
pt.data$lr.psad.miss<-as.numeric(is.na(demo.data$PSADatDX))
#table(demo.data$PSAD_criteria)
#xtabs(~lr.psad+lr.psad.miss, data=pt.data)

pt.data$lr.psad[pt.data$lr.psad.miss==1]<-1

# numeric DOB
#### This takes dates with a 2-digit year of for DOB
d <- as.Date(demo.data$DOB, "%m/%d/%y")
pt.data$dob.num <- as.numeric(as.Date(ifelse(d > "2000-1-1", format(d, "19%y-%m-%d"), format(d))))
data.check(condition=as.logical(sum(is.na(pt.data$dob.num))==0), message="Problems converting patients' date of birth. NAs occurred.")
#### I haven't included any data checks that DOBs are reasonable


# numeric date of death
pt.data$dod.num <- as.numeric(as.Date(demo.data$Date_of_death, "%m/%d/%y"))
# will be NA if patient did not die
data.check(condition=as.logical(sum(demo.data$Date_of_death=="")==sum(is.na(pt.data$dod.num))), message="Problems converting patients' date of death. NAs occurred for patients with non-missing DOD.")
#### I haven't included any data checks that DODs are reasonable



### 3. Identify diagnostic biopsy for all patients. Limit analysis to patients with GS<7 at dx. Put bx dx date (and age) into pt.data.------------
#define vector used to indicate records for removal
pt.data$rm <- rep(0,n)

n_bx <- dim(bx.data)[1]
bx.data$rm <- rep(0, n_bx)

#remove patients without diagnostic biopsy.
#### This code does not check whether patients without diagnostic biopsies have later biopsies. (Doesn't seem to be the case in jan 2016 data pull.) All patients without dx bx are excluded from the analysis.
rm.id <- pt.data$id[! pt.data$id %in% bx.data$id[bx.data$Which_Biopsy=="D"]]
if(length(rm.id)>0){
	warning(paste0(length(rm.id), " patients were removed from the analysis because they did not have a diagnostic biopsy."))}

pt.data$rm[pt.data$id%in%rm.id] <- 1
pt.data <- pt.data[pt.data$rm==0,]

bx.data$rm[bx.data$id %in% rm.id] <- 1
bx.data <- bx.data[bx.data$rm==0,]
n_bx <- dim(bx.data)[1]

#remove pts with dx bx GS>6
#### In my version of the data, missing GS are coded with "" not NA
bx.data$gleason[is.na(bx.data$gleason)] <- "" # 110622 Zitong recode missing gleason to keep consistency with Coley
warning(paste0("Gleason score is missing for ", sum(bx.data$gleason==""), " biopsies. ", sum(bx.data$NumPosCores>0 & bx.data$gleason=="" & !is.na(bx.data$NumPosCores)), " of these biopsies had at least one positive core found. Analysis continues under the assumption that none of these biopsies were grade reclassifications."))

#make a variable for grade reclassification
bx.data$rc <- rep(1,n_bx)
bx.data$rc[bx.data$gleason==6] <- 0
bx.data$rc[bx.data$gleason==""] <- NA

#remove those with grade RC at dx bx
rm.id <- bx.data$id[bx.data$rc==1 & !is.na(bx.data$rc) & bx.data$Which_Biopsy=="D"]
warning(paste0("Gleason score at diagnostic biopsy is > 6 for ", length(rm.id), " patients. These patients will be removed from the analysis and no PPGG predictions will be available for them. "))

pt.data$rm[pt.data$id %in% rm.id] <- 1
pt.data <- pt.data[pt.data$rm==0,]
bx.data$rm[bx.data$id %in% rm.id] <- 1
bx.data <- bx.data[bx.data$rm==0,]

#update n
n <- dim(pt.data)[1]
n_bx <- dim(bx.data)[1]


### Remove patients with fewer than 2 biopsies (i.e., only a diagnostic biopsy)
pt.data$num.bx<-rep(0,n)
for(i in 1:n){
  pt.data$num.bx[i]<-sum(
    bx.data$clinical_PTnum==pt.data$clinical_PTnum[i])  }
table(pt.data$num.bx)

pt.data$rm[pt.data$num.bx<2]<-1
bx.data$rm[bx.data$id %in% pt.data$id[pt.data$rm==1]] <- 1
pt.data <- pt.data[pt.data$rm==0,]
bx.data <- bx.data[bx.data$rm==0,]
#update n
n <- dim(pt.data)[1]
n_bx <- dim(bx.data)[1]


#add volume criteria at dx
pt.data$lr.vol<-pt.data$lr.vol.miss<-rep(0,n)
for(i in 1:n){
  if(is.na(bx.data$NumPosCores[bx.data$id==pt.data$id[i] & bx.data$Which_Biopsy=="D"])
      | is.na(bx.data$MaxPrcntCore[bx.data$id==pt.data$id[i] & bx.data$Which_Biopsy=="D"])){
    pt.data$lr.vol.miss[i]<-1}
  else{ if(bx.data$NumPosCores[bx.data$id==pt.data$id[i] & bx.data$Which_Biopsy=="D"]>2
           | bx.data$MaxPrcntCore[bx.data$id==pt.data$id[i] & bx.data$Which_Biopsy=="D"]>50){
    pt.data$lr.vol[i]<-1
  }}
}

#table(pt.data$lr.vol)
#table(pt.data$lr.vol.miss)

#patients with missing volume information are considered low risk (high volume) vs. very low risk
pt.data$lr.vol[pt.data$lr.vol.miss==1]<-1

#xtabs(~vlr+lr.vol, data=pt.data)
#xtabs(~vlr+lr.vol.miss, data=pt.data)




#numeric date for all bx; numeric date of dx bx into pt.data
head(bx.data$Biopsy_Date)
# bx.data$bx.dt.num <- as.numeric(as.Date(bx.data$Biopsy_Date,"%d-%b-%y"))
bx.data$bx.dt.num <- as.numeric(as.Date(bx.data$Biopsy_Date,"%m/%d/%y")) #112522 Zitong fix biopsy date format for current dataset
if(sum(is.na(bx.data$bx.dt.num))>0){
  stop("Incorrect biopsy data Biopsy_Date format")}
pt.data$dx.dt.num <- rep(0,n)
for(i in 1:n){
  pt.data$dx.dt.num[i] <- bx.data$bx.dt.num[bx.data$id==pt.data$id[i] & bx.data$Which_Biopsy=="D"]}
##### could put some checks here- e.g. bx can't occur after death

#patient age at dx
pt.data$dx.age <- (pt.data$dx.dt.num-pt.data$dob.num)/365




### 4. Treatment data: Remove pts removed from pt.data; remove records without date; only keep 1st tx per pt. Put tx date, type of tx, and post-surgery GS in pt.data ----------------
#remove pts removed from pt.data
tx.data <- tx.data[tx.data$id %in% pt.data$id,]

#remove records without date
#### here, missing treatment date is " ", not NA
if(sum(tx.data$Type_Of_Treatment=="Surg" & tx.data$Treatment_Date=="")>0){
	warning(paste0(sum(tx.data$Type_Of_Treatment=="Surg" & tx.data$Treatment_Date==""), " patients do not a have date for surgery received. As a result, that surgery was performed for these patients can not be included in the analysis."))}
##### could also include warning message for other patients with no tx date

tx.data <- tx.data[!tx.data$Treatment_Date=="",]
n_tx <- dim(tx.data)[1]

# make numeric date for treatments
tx.data$tx.dt.num <- as.numeric(as.Date(tx.data$Treatment_Date, "%m/%d/%y"))

#put first tx date into pt.data; remove later dates
pt.data$tx.dt.num <- rep(NA,n)
tx.data$rm <- rep(0,n_tx)
for(i in 1:n){
	if(pt.data$id[i] %in% tx.data$id){
		pt.data$tx.dt.num[i] <- min(tx.data$tx.dt.num[tx.data$id==pt.data$id[i]])
		tx.data$rm[tx.data$id==pt.data$id[i] & !tx.data$tx.dt.num==pt.data$tx.dt.num[i]] <- 1
		}}
tx.data <- tx.data[tx.data$rm==0,]
n_tx <- dim(tx.data)[1]

#calculate time-until-treatment for patients
pt.data$time.until.tx <- (pt.data$tx.dt.num - pt.data$dx.dt.num)/365
data.check(condition=as.logical(sum(!is.na(pt.data$time.until.tx) & pt.data$time.until.tx<0)==0), message="Some patients have initial treatment dates occurring prior to diagnostic biopsy.")

#indicate patients with surgery in pt.data
tx.data$Type_Of_Treatment[tx.data$Type_Of_Treatment=="rrp"]<-"Surg"
pt.data$surg <- rep(0,n)
for(i in 1:n){
	if(pt.data$id[i] %in% tx.data$id[tx.data$Type_Of_Treatment=="Surg"]){
		pt.data$surg[i] <- 1 }	}


# group post-surgery pathological GS into prognostic grade group
tx.data$true.pgg<-rep(NA,n_tx)
tx.data$true.pgg[tx.data$DominantRRP_Grade=="3+2" | tx.data$DominantRRP_Grade=="2+3" | tx.data$DominantRRP_Grade=="3+3"]<-1
tx.data$true.pgg[tx.data$DominantRRP_Grade=="3+4"]<-2
tx.data$true.pgg[tx.data$DominantRRP_Grade=="4+3"]<-3
tx.data$true.pgg[tx.data$DominantRRP_Grade=="4+4" | tx.data$DominantRRP_Grade=="4+5" | tx.data$DominantRRP_Grade=="5+4" | tx.data$DominantRRP_Grade=="5+5"]<-4

if(sum(tx.data$DominantRRP_Grade=="" & tx.data$Type_Of_Treatment=="Surg")>0){
	warning(paste0(sum(tx.data$DominantRRP_Grade=="" & tx.data$Type_Of_Treatment=="Surg"), " patients had surgery performed but did not have a post-surgery pathological grade determination."))}

# record true PGG in pt.data (and binary PGG, PGG=1 or PGG>1)
pt.data$true.pgg<-pt.data$true.pgg.bin<-rep(NA,n)
for(i in 1:n){
	if(pt.data$surg[i]==1){
		pt.data$true.pgg[i] <- tx.data$true.pgg[tx.data$id==pt.data$id[i]]
		if(!is.na(pt.data$true.pgg[i])){
			pt.data$true.pgg.bin[i] <- as.numeric(pt.data$true.pgg[i]>1)} } }



### 5. PSA data- remove patients removed from pt.data; remove records without value or date; make PSA date numeric; calculate time since diagnosis and age at PSA; remove PSA dates after treatment; remove patients with less than 2 PSA tests; log-transform PSA; pull prostate volume data from bx.data; record date of last PSA test ----------
#remove PSA observations more than 1 year prior to dx bx
psa.data <- psa.data[psa.data$id %in% pt.data$id,]
psa.data <- psa.data[!is.na(psa.data$Total_PSA),]
psa.data <- psa.data[!is.na(psa.data$PSA_Date),]
n_psa <- dim(psa.data)[1]

#numeric PSA dates
psa.data$dt.num <- as.numeric(as.Date(psa.data$PSA_Date, "%m/%d/%y"))
psa.data <- psa.data[!is.na(psa.data$dt.num),] ## 080821
n_psa <- dim(psa.data)[1]## 080821
data.check(condition=as.logical(sum(is.na(psa.data$dt.num))==0), message="Problems converting PSA test date. NAs occurred.")

#time since dx, age for each PSA test
psa.data$time.since.dx <- psa.data$age <- vector(length=n_psa)
for(i in 1:n){
	psa.data$time.since.dx[psa.data$id==pt.data$id[i]] <- (psa.data$dt.num[psa.data$id==pt.data$id[i]] - pt.data$dx.dt.num[i])/365
	psa.data$age[psa.data$id==pt.data$id[i]] <- (psa.data$dt.num[psa.data$id==pt.data$id[i]] - pt.data$dob[i])/365}

data.check(condition=as.logical(sum(psa.data$age==0)==0), message="Some PSA tests occurred prior to DOB.")


#remove PSA dates after treatment
psa.data$rm <- rep(0,n_psa)
for(i in 1:n){
	if(!is.na(pt.data$time.until.tx[i])){
		psa.data$rm[psa.data$id==pt.data$id[i]] <- as.numeric(psa.data$time.since.dx[psa.data$id==pt.data$id[i]] > pt.data$time.until.tx[i])} }

psa.data <- psa.data[psa.data$rm==0,]
n_psa <- dim(psa.data)[1]
table(psa.data$rm) #no PSA after treatment

#remove PSA tests more than 1 year prior to diagnosis
sum(psa.data$time.since.dx< (-1) )
mean(psa.data$time.since.dx< (-1) ) #0.20

length(unique(psa.data$id[psa.data$time.since.dx<(-1)])) #
length(unique(psa.data$id[psa.data$time.since.dx<(-1)]))/n #

psa.data<-psa.data[psa.data$time.since.dx>= (-1),]


#remove patients with less than 2 PSA tests
for(i in 1:n){
	if(sum(psa.data$id==pt.data$id[i])<2){
		pt.data$rm[i] <- 1
		psa.data$rm[psa.data$id==pt.data$id[i]] <- 1 } }
table(pt.data$rm)
pt.data <- pt.data[pt.data$rm==0,]
n <- dim(pt.data)[1]
psa.data <- psa.data[psa.data$rm==0,]
n_psa <- dim(psa.data)[1]


# log-transform PSA
psa.data$log.psa <- log(psa.data$Total_PSA)# + 0.1)


# pull prostate volume data from bx.data
# use average prostate volume from all biopsies
pt.data$avg.vol <- vector(length=n)
for(i in 1:n){
	pt.data$avg.vol[i] <- mean(bx.data$VolatBX[bx.data$id==pt.data$id[i]], na.rm=T)}
if(sum(is.na(pt.data$avg.vol))>0){
	warning(paste0( sum(is.na(pt.data$avg.vol)), " patients did not have prostate volume assessment at any biopsy. Average prostate volume will be assumed for these patients."))}

# mean- variance- standardized average prostate volume
pt.data$std.vol[!is.na(pt.data$avg.vol)] <- scale(pt.data$avg.vol[!is.na(pt.data$avg.vol)])
pt.data$std.vol[is.na(pt.data$std.vol)] <- 0

#add this to PSA data
psa.data$std.vol <- vector(length=n_psa)
for(i in 1:n){
	psa.data$std.vol[psa.data$id==pt.data$id[i]] <- pt.data$std.vol[i] }

# record date of last PSA test
pt.data$last.psa.dt.num <- vector(length=n)
for(i in 1:n){
	pt.data$last.psa.dt.num[i] <- max(psa.data$dt.num[psa.data$id==pt.data$id[i]])}




### 6 Biopsy data- remove records from patients removed due to PSA criteria; Remove any records with missing biopsy date; Calculate time since diagnosis and remove observations prior to dx or after tx; Define bx PGG results; Move RC data, RC and final biopsy date to pt.data; remove biopsies that occur after initial reclassification; tidy NPC and MPC data ------------

#remove records from patients removed due to PSA criteria
bx.data <- bx.data[bx.data$id%in%pt.data$id,]
n_bx <- dim(bx.data)[1]

#remove records with missing date
bx.data <- bx.data[!is.na(bx.data$bx.dt.num),]
n_bx <- dim(bx.data)[1]

#add avg, std prostate volume
bx.data$avg.vol<-bx.data$std.vol <- vector(length=n_bx)
for(i in 1:n){
  bx.data$std.vol[bx.data$id==pt.data$id[i]] <- pt.data$std.vol[i]
  bx.data$avg.vol[bx.data$id==pt.data$id[i]] <- pt.data$avg.vol[i] }
#summary(bx.data$std.vol)
#summary(bx.data$avg.vol)
#hist(bx.data$avg.vol, breaks=20)


#calculate time since dx for each bx
bx.data$time.since.dx <- vector(length=n_bx)
for(i in 1:n){
	bx.data$time.since.dx[bx.data$id==pt.data$id[i]] <- (bx.data$bx.dt.num[bx.data$id==pt.data$id[i]] - pt.data$dx.dt.num[i])/365}

if(sum(bx.data$time.since.dx<0)>0){
	warning(paste0( length(unique(bx.data$id[bx.data$time.since.dx<0])), " patients had biopsies performed prior to diagnosis. These biopsy records will be removed from the analysis and analysis will proceed. Check for data entry error."))}

bx.data <- bx.data[bx.data$time.since.dx>=0,]
n_bx <- dim(bx.data)[1]

#remove any bx after tx
for(i in 1:n){
	if(!is.na(pt.data$time.until.tx[i])){
		bx.data$rm[bx.data$id==pt.data$id[i]] <- as.numeric(bx.data$time.since.dx[bx.data$id==pt.data$id[i]] > pt.data$time.until.tx[i])}}

if(sum(bx.data$rm)>0){
	print()
	warning(paste0(length(unique(bx.data$id[bx.data$rm==1])), " patients had biopsies performed after treatment. These biopsy records will be removed from the analysis and analysis will proceed. Check for data entry error."))}

bx.data <- bx.data[bx.data$rm==0,]
n_bx <- dim(bx.data)[1]
#### Current script doesn't check that biopsies don't occur after death




#biopsy results- PGG
#check biopsy gleason vs. BXshowCa
if(sum(bx.data$BXshowCa==1 & bx.data$gleason=="")>0){
  warning(paste0(sum(bx.data$BXshowCa==1 & bx.data$gleason==""), " records have BXshowCa==1 but missing biopsy Gleason score. These records cannot be included in the Biopsy PGG analysis since this information is missing."))}

#create biopsy- and patient-level biopsy PGG variable
bx.data$pgg <- rep(NA,n_bx)
bx.data$pgg[bx.data$gleason=="6"] <- 1
bx.data$pgg[bx.data$gleason=="7(3+4)"] <- 2
bx.data$pgg[bx.data$gleason=="7(4+3)"] <- 3
bx.data$pgg[bx.data$gleason=="8(4+4)" | bx.data$gleason=="8(3+5)"
            | bx.data$gleason=="8(5+3)" | bx.data$gleason=="9(4+5)"
            | bx.data$gleason=="9(5+4)" | bx.data$gleason=="10(5+5)"] <- 4

#move RC, PGG data to pt.data
pt.data$rc <- rep(0,n)
for(i in 1:n){
	if(max(bx.data$rc[bx.data$id==pt.data$id[i]], na.rm=T)==1){pt.data$rc[i]<-1} }

pt.data$bx.pgg <- rep(0,n)
for(i in 1:n){
	pt.data$bx.pgg[i] <- max(bx.data$pgg[bx.data$id==pt.data$id[i]], na.rm=T)}



#date of first rc
pt.data$rc.dt.num <- rep(NA,n)
for(i in 1:n){
	if(pt.data$rc[i]==1){
		pt.data$rc.dt.num[i] <- min(bx.data$bx.dt.num[
		  bx.data$id==pt.data$id[i] & bx.data$rc==1 & !is.na(bx.data$rc)]) } }




#date of last bx
pt.data$last.bx.dt.num <- vector(length=n)
for(i in 1:n){
	pt.data$last.bx.dt.num[i] <- max(bx.data$bx.dt.num[bx.data$id==pt.data$id[i]])}

#remove biopsies that occur after initial reclassification
#### Since only 17 patients had additional biopsies after initial reclassification (and only 3 had surgery), I exclude additional biopsies performed after initial grade reclassification. I hope to life this exclusion when we have data on more patients who continue with biopsies afer upgrading.
for(i in 1:n){
	if(pt.data$rc[i]==1){
		if(!pt.data$rc.dt.num[i]==pt.data$last.bx.dt.num[i]){
			bx.data$rm[bx.data$id==pt.data$id[i] & bx.data$bx.dt.num>pt.data$rc.dt.num[i]] <- 1
			pt.data$last.bx.dt.num[i] <- pt.data$rc.dt.num[i] } } }

if(sum(bx.data$rm)>0){
  warning(paste0(length(unique(bx.data$id[bx.data$rm==1])), " patients have a total of ", sum(bx.data$rm==1), " biopsies after initial reclassification. These biopsies are removed from this analysis. In the future, we plan to extend the model to accomodate biopsies after initial reclassification."  ))}

bx.data <- bx.data[bx.data$rm==0,]
n_bx <- dim(bx.data)[1]

#time until RC
pt.data$time.until.rc <- vector(length=n)
for(i in 1:n){
	if(pt.data$rc[i]==1){
	  ## warning occurs for i = 1247, set RHS to be unique values _ Zitong
		pt.data$time.until.rc[i] <-
		  unique(bx.data$time.since.dx[bx.data$id==pt.data$id[i] & bx.data$rc==1 & !is.na(bx.data$rc)])	}
		else{pt.data$time.until.rc[i]<-NA}}




###NPC and MPC data

##rename NPC and MPC
#rename variable "npc" for brevity and to save original data
bx.data$npc <- bx.data$NumPosCore

#rename variable "MPC" for brevity and to save original data
bx.data$mpc <- bx.data$MaxPrcntCore


##examine diagnosis values
#deal with missing values at diagnosis
if(sum(is.na(bx.data$npc) & bx.data$Which_Biopsy=="D")>0){
	warning(paste0(sum(is.na(bx.data$npc) & bx.data$Which_Biopsy=="D"), " patients have missing # of positive cores for diagnostic biopsy. One positive core at diagnosis is assumed. This information is only used to predict future decisions to have biopsy or surgery. (This patient is separately categorized as Low Risk- Volume for having missing volume data at diagnosis.)"))}
bx.data$npc[is.na(bx.data$npc) & bx.data$Which_Biopsy=="D"] <- 1

#npc should be at least 1 at dx bx
if(sum(bx.data$npc==0 & !is.na(bx.data$npc) & bx.data$Which_Biopsy=="D")>0){
	warning(paste0(sum(bx.data$npc==0 & !is.na(bx.data$npc) & bx.data$Which_Biopsy=="D"), " patients have zero positive cores at diagnostic biopsy. One positive core at diagnosis is assumed. This information is only used to predict future decisions to have biopsy or surgery. Check for data entry error"))}
bx.data$npc[bx.data$npc==0 & !is.na(bx.data$npc) & bx.data$Which_Biopsy=="D"] <- 1

#deal with missing values
if(sum(is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D")>0){
  warning(paste0(sum(is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D"), " patients have missing maximum % positive for diagnostic biopsy. 1% positive at diagnosis is assumed. This information is only used to predict future decisions to have biopsy or surgery. (This patient is separately categorized as Low Risk- Volume for having missing volume data at diagnosis.)"))}
bx.data$mpc[is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D"] <- 1

#mpc should be at least 1 at dx bx
if(sum(bx.data$mpc==0 & !is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D")>0){
  warning(paste0(sum(bx.data$mpc==0 & !is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D"), " patient ids above have 0 maximum % positive at diagnostic biopsy. 1% positive at diagnosis is assumed. This information is only used to predict future decisions to have biopsy or surgery. Check for data entry error."))}
bx.data$npc[bx.data$mpc==0 & !is.na(bx.data$mpc) & bx.data$Which_Biopsy=="D"] <- 1




#NPC uses the BXshowCa variable, so we want to check accuracy here
table(bx.data$gleason[bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa)])
if(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) #& #!bx.data$gleason=="6"
        & !bx.data$gleason=="")>0){
  warning(paste0(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) #&  !bx.data$gleason=="6"
                       & !bx.data$gleason==""),
                 " records have BXshowCa==0 but Gleason score >= 6. These records will be changed to BXshowCa==1. Check for data entry error"))}
if(sum(is.na(bx.data$BXshowCa) # & !bx.data$gleason=="6"
       & !bx.data$gleason=="")>0){
  warning(paste0(sum(is.na(bx.data$BXshowCa) #&!bx.data$gleason=="6"
                       & !bx.data$gleason==""),
                 " records have is.na(BXshowCa)==1 but Gleason score >= 6. These records will be changed to BXshowCa==1. Check for data entry error"))}

table(bx.data$npc[bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa)])
if(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
       bx.data$npc>0 & !is.na(bx.data$npc))>0){
  warning(paste0(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                       bx.data$npc>0 & !is.na(bx.data$npc)),
                 " records have BXshowCa==0 but NumPosCore > 0. These records will be changed to BXshowCa==1. Check for data entry error"))}
if(sum(is.na(bx.data$BXshowCa) &
       bx.data$npc>0 & !is.na(bx.data$npc))>0){
  warning(paste0(sum(is.na(bx.data$BXshowCa) &
                       bx.data$npc>0 & !is.na(bx.data$npc)),
                 " records have is.na(BXshowCa)==1 but NumPosCore > 0. These records will be changed to BXshowCa==1. Check for data entry error"))}


table(bx.data$mpc[bx.data$BXshowCa==0])
if(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
       bx.data$mpc>0 & !is.na(bx.data$mpc))>0){
  warning(paste0(sum(bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                       bx.data$mpc>0 & !is.na(bx.data$mpc)),
                 " records have BXshowCa==0 but MaxPrcntCore > 0. These records will be changed to BXshowCa==1. Check for data entry error"))}
if(sum(is.na(bx.data$BXshowCa) &
       bx.data$mpc>0 & !is.na(bx.data$mpc))>0){
  warning(paste0(sum(is.na(bx.data$BXshowCa) &
                       bx.data$mpc>0 & !is.na(bx.data$mpc)),
                 " records have is.na(BXshowCa)==1 but MaxPrcntCore > 0. These records will be changed to BXshowCa==1. Check for data entry error"))}



#update records
table(bx.data$BXshowCa)
bx.data$BXshowCa[bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) #& !bx.data$gleason=="6"
                    & !bx.data$gleason==""]<-1
bx.data$BXshowCa[is.na(bx.data$BXshowCa) # & !bx.data$gleason=="6"
                    & !bx.data$gleason==""]<-1
bx.data$BXshowCa[bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                   bx.data$npc>0 & !is.na(bx.data$npc)]<-1
bx.data$BXshowCa[is.na(bx.data$BXshowCa) &
                   bx.data$npc>0 & !is.na(bx.data$npc)]<-1
bx.data$BXshowCa[bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                   bx.data$mpc>0 & !is.na(bx.data$mpc)]<-1
bx.data$BXshowCa[is.na(bx.data$BXshowCa) &
                   bx.data$mpc>0 & !is.na(bx.data$mpc)]<-1
table(bx.data$BXshowCa)


#is BXshowCa==0 used correctly?
if(sum(bx.data$npc>0 & !is.na(bx.data$npc) &
       bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
       !bx.data$Which_Biopsy=="D")>0){
  warning(paste0(sum(bx.data$npc>0 & !is.na(bx.data$npc) &
                       bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                       !bx.data$Which_Biopsy=="D"),
                 " records have >0 positive cores but are marked as BXshowCa==0. NumPosCores value will be used. Check for data entry error."))}

if(sum(bx.data$mpc>0 & !is.na(bx.data$mpc) &
       bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
       !bx.data$Which_Biopsy=="D")>0){
  warning(paste0(sum(bx.data$mpc>0 & !is.na(bx.data$mpc) &
                       bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
                       !bx.data$Which_Biopsy=="D"),
                 " records have >0 max % involvement but are marked as BXshowCa==0. MaxPrcntCore value will be used. Check for data entry error."))}



#tidy NPC post-diagnosis
#table(bx.data$npc[!bx.data$Which_Biopsy=="D"]) #0-12 positive cores
#sum(is.na(bx.data$npc) & !bx.data$Which_Biopsy=="D") #27 with missing NPC
#table(bx.data$BXshowCa[is.na(bx.data$npc) & !bx.data$Which_Biopsy=="D"])  #12 missing NPC don't show cancer
#sum(is.na(bx.data$BXshowCa)) #some missingness here as well


if(sum(is.na(bx.data$npc) & bx.data$BXshowCa==0)>0){
  warning(paste0(sum(is.na(bx.data$npc) & bx.data$BXshowCa==0),
                 " records have missing data for NumPosCores and BXshowCa==0. These records will be assigned NumPosCores==0"))}

bx.data$npc[is.na(bx.data$npc) &
              bx.data$BXshowCa==0 & !is.na(bx.data$BXshowCa) &
              !bx.data$Which_Biopsy=="D"]<-0

if(sum(is.na(bx.data$npc) & bx.data$BXshowCa==1 & !is.na(bx.data$BXshowCa))>0){
  warning(paste0(sum(is.na(bx.data$npc) & bx.data$BXshowCa==1),
                 " records have missing data for NumPosCores and BXshowCa==1. These records will be excluded from the # positive cores regression component of the prediction model. These records will be included in the biopsy PGG regression component if biopsy grades are available."))}

if(sum(bx.data$npc==0 & !is.na(bx.data$npc)
       & !bx.data$gleason=="")>0){
  warning(paste0(sum(bx.data$npc==0 & !is.na(bx.data$npc)
                     & !bx.data$gleason==""),
                 " records have NumPosCores==0 but a Gleason score >=6. These records will be excluded from the # positive cores regression component of the prediction model. Check for data entry error."))}
bx.data$npc[bx.data$npc==0 & !is.na(bx.data$npc) & !bx.data$gleason==""] <- NA
#xtabs(~bx.data$npc + bx.data$gleason)



#MPC
#table(bx.data$mpc[!bx.data$Which_Biopsy=="D"]) #1-100% positive cores
#sum(is.na(bx.data$mpc) & !bx.data$Which_Biopsy=="D") #2139 with missing NPC
#table(bx.data$BXshowCa[is.na(bx.data$mpc) & !bx.data$Which_Biopsy=="D"])  #12 missing NPC don't show cancer
#table(bx.data$npc[is.na(bx.data$mpc) & !bx.data$Which_Biopsy=="D"])  #12 missing NPC don't show cancer

if(sum(is.na(bx.data$mpc) & bx.data$npc==0 & !is.na(bx.data$npc))>0){
  warning(paste0(sum(is.na(bx.data$mpc) & bx.data$npc==0 & !is.na(bx.data$npc)),
                 " records have missing data for MaxPrcntCore and NumPosCores==0. These records will be assigned MaxPrcntCore==0"))}

bx.data$mpc[is.na(bx.data$mpc) & !is.na(bx.data$npc) & bx.data$npc==0 &
              !bx.data$Which_Biopsy=="D"]<-0

#sum(is.na(bx.data$mpc))
#table(bx.data$npc[is.na(bx.data$mpc)])
#sum(is.na(bx.data$mpc) & is.na(bx.data$npc))

if( sum(is.na(bx.data$mpc) & bx.data$npc>0 & !is.na(bx.data$npc)) >0){
  warning(paste0(sum(is.na(bx.data$mpc) & bx.data$npc>0 & !is.na(bx.data$npc) ),
                 " records have missing data for MaxPrcntCore and NumPosCores>0. These records will be excluded from the max % involvement regression component of the prediction model."))}

if(sum(is.na(bx.data$mpc) & is.na(bx.data$npc))>0){
  warning(paste0(sum(is.na(bx.data$mpc) & is.na(bx.data$npc) ),
                 " records have missing data for MaxPrcntCore and NumPosCores. These records will be excluded from the max % involvement regression component of the prediction model. "))}



##MRI-targeted biopsy
#indicate bx where MRI was performed

#table(bx.data$mritargetedbx) #values 0-3
#sum(is.na(bx.data$mritargetedbx)) #10
bx.data$mri<-rep(0,n_bx)
bx.data$mri[bx.data$mritargetedbx>0 & !is.na(bx.data$mritargetedbx)]<-1
#table(bx.data$mri)
if(sum(is.na(bx.data$mritargetedbx))>0){
  warning(paste0(sum(is.na(bx.data$mritargetedbx)),
                 " biopsies are missing information for mritargetedbx. We will assume that the biopsies were not performed with MRI.")) }

#MRI-targeted biopsy assign specific pirads score to mri targeted bxs 1126 Zitong
if(sum(bx.data$max_pirads>5, na.rm = T)>0){
  warning(paste0(sum(bx.data$max_pirads>5, na.rm = T),
                 " biopsies have max pirads > 5. Check data entry error. ")) }



##Number of cores sampled
#rename for brevity and so that original data are preserved
bx.data$ncs<-bx.data$NumCoreSampled
#table(bx.data$ncs)
#hist(bx.data$ncs)

#assign likely number of cores where data are missing
#bx.data$ncs[bx.data$ncs == 0] <- NA #110622 Zitong - recode NA values to ncs = 0, for adapting new data
sum(is.na(bx.data$ncs))
ncs6.date<-as.numeric(as.Date("1/1/03","%m/%d/%y"))

if(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num < ncs6.date)>0){
  warning(paste0(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num < ncs6.date),
                 " records have missing NumCoreSampled and were performed prior to 1/1/2012. We assume 6 cores sampled."))}

if(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num >= ncs6.date & bx.data$mri==0)>0){
  warning(paste0(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num >= ncs6.date & bx.data$mri==0),
                 " records have missing NumCoreSampled, were performed on or after 1/1/2012, and were not MRI-targeted. We assume 12 cores sampled."))}


if(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num >= ncs6.date & bx.data$mri==1)>0){
  warning(paste0(sum(is.na(bx.data$ncs) & bx.data$bx.dt.num >= ncs6.date & bx.data$mri==1),
                 " records have missing NumCoreSampled, were performed on or after 1/1/2012, and were MRI-targeted. We assume 15 cores sampled."))}

bx.data$ncs[is.na(bx.data$ncs)& bx.data$bx.dt.num < ncs6.date] <- 6
bx.data$ncs[is.na(bx.data$ncs)& bx.data$bx.dt.num >= ncs6.date
            & bx.data$mri==0] <- 12
bx.data$ncs[is.na(bx.data$ncs)& bx.data$bx.dt.num >= ncs6.date
            & bx.data$mri==1] <- 15


### laterality
#has cancer been detected on both sides of the prostate?
#table(bx.data$lateralitybx)
bx.data$lat<-bx.data$lateralitybx

#xtabs(~BXshowCa+lat,data=bx.data) #are these cumulative laterality and bx after bilateral disease was found?
#xtabs(~NumPosCores+lat,data=bx.data)
#xtabs(~npc+lat,data=bx.data)
#table(bx.data$lat[is.na(bx.data$NumPosCores) & !bx.data$Which_Biopsy=="D"])

pt.data$lat<-rep(0,n)
for(i in 1:n){
  if(max(bx.data$lat[bx.data$id==pt.data$id[i]])==1){
    pt.data$lat[i]<-1}}
#table(pt.data$lat) #692 with bilateral pca

for(i in 1:n){
  if(pt.data$lat[i]==1){
    lat.time<-min(bx.data$time.since.dx[bx.data$lat==1 & bx.data$id==pt.data$id[i]])
    bx.data$lat[bx.data$id==pt.data$id[i] & bx.data$time.since.dx>lat.time]<-NA
    lat.time<-NULL } }

#table(bx.data$lat)
#sum(is.na(bx.data$lat))
#xtabs(~BXshowCa+lat,data=bx.data)








### 7. Calculate total follow-up time. Categorize patients based on outcomes, LTF ----------

pt.data$time.fup <- rep(0,n) #total time of follow-up (use for making annual interval dataframe)
pt.data$active <- rep(0,n) # have not had RC or surgery; have had bx or PSA within 2 years of data pull
pt.data$ltf <- rep(0,n) # have not had RC or surgery; lost to follow-up, e.g. more than 2 years between last bx/psa and data pull
pt.data$rc.censored <- rep(0,n) # RC observed but no surgery up until time of other treatment or 2 years post-RC follow-up
pt.data$rc.active <- rep(0,n) # RC observed less than 2 years before data pull, patient has not yet elected treatment

#max time until surgery after rc
max.time.until.surg <- ceiling(max(pt.data$time.until.tx[!is.na(pt.data$time.until.tx) & !is.na(pt.data$time.until.rc)] - pt.data$time.until.rc[!is.na(pt.data$time.until.tx) & !is.na(pt.data$time.until.rc)]))

for(i in 1:n){
#follow-up time for those with treatment
	if(!is.na(pt.data$tx.dt.num[i])){
		pt.data$time.fup[i] <- (pt.data$tx.dt.num[i] - pt.data$dx.dt.num[i]) / 365}
	else{
#follow-up time for those with reclassification but no treatment yet
		if(pt.data$rc[i]==1){
#maximum follow-up time is 2 years since rc, compare to date of data pull or death
			max.fup <- pt.data$rc.dt.num[i] + (365*max.time.until.surg)
			if(is.na(pt.data$dod.num[i])){ #for those that did not die
				if(max.fup > date.pull){
					pt.data$rc.active[i] <- 1
					pt.data$time.fup[i] <- (date.pull - pt.data$dx.dt.num[i]) / 365}
				if(max.fup <= date.pull){
					pt.data$rc.censored[i] <- 1
					pt.data$time.fup[i] <- (max.fup - pt.data$dx.dt.num[i]) / 365} }
					else{ #died after RC, prior to treatment
				if(max.fup > pt.data$dod.num[i]){
					pt.data$time.fup[i] <- (pt.data$dod.num[i] - pt.data$dx.dt.num[i]) / 365}
				if(max.fup <= pt.data$dod.num[i]){
					pt.data$rc.censored[i] <- 1
					pt.data$time.fup[i] <- (max.fup - pt.data$dx.dt.num[i]) / 365} }
		}
		else{ #follow-up time for those without reclassification or treatment
			max.fup <- max(pt.data$last.bx.dt.num[i], pt.data$last.psa.dt.num[i]) + (365*2)
			if(is.na(pt.data$dod.num[i])){ #pts prior to RC and death
				if(max.fup > date.pull){
					pt.data$active[i] <- 1
					pt.data$time.fup[i] <- (date.pull - pt.data$dx.dt.num[i]) / 365}
				if(max.fup <= date.pull){
					pt.data$ltf[i] <- 1
					pt.data$time.fup[i] <- (max.fup - pt.data$dx.dt.num[i]) / 365}
					}
					else{ #pts who died prior to RC
				if(max.fup > pt.data$dod.num[i]){
					pt.data$time.fup[i] <- (pt.data$dod.num[i] - pt.data$dx.dt.num[i]) / 365}
				if(max.fup <= pt.data$dod.num[i]){
					pt.data$ltf[i] <- 1
					pt.data$time.fup[i] <- (max.fup - pt.data$dx.dt.num[i]) / 365} }
					}		}
					}

data.check(condition=as.logical(sum(is.na(pt.data$time.fup))==0), message="There was a problem identifying total follow-up time for one or more patients. Email Yates; she will check original script.")
data.check(condition=as.logical(sum(!is.na(pt.data$time.until.tx) & !pt.data$time.until.tx==pt.data$time.fup)==0), message="Follow-up time does not match time-until-tx for at least one patient. Email Yates; she will check orginal script.")





### 8. Order patients based on observed eta --------
#### The need to order patient data in this way is an artifact of JAGS indexing requirements
ordered<-order(pt.data$true.pgg)
pt.data<-pt.data[ordered,]
pt.data$subj<-c(1:n)

psa.data$subj<-rep(0,n_psa)
bx.data$subj<-rep(0,n_bx)
for(i in 1:n){
	psa.data$subj[psa.data$id==pt.data$id[i]] <- i
	bx.data$subj[bx.data$id==pt.data$id[i]] <- i }



### 9. Shaping biopsy data. Put data into annual intervals ----------
#112622 Zitong - this section is modified to include mritargeted bx pirads score for new data

#create dataframe with one row per patient per year under observation, including dx bx
bx.subj <- rep(1, (ceiling(pt.data$time.fup[pt.data$subj==1]) + 1) )
bx.time <- c(0:ceiling(pt.data$time.fup[pt.data$subj==1]))
for(i in 2:n){
	bx.subj <- c(bx.subj, rep(i, (ceiling(pt.data$time.fup[pt.data$subj==i]) + 1)) )
	bx.time <- c(bx.time, c(0:ceiling(pt.data$time.fup[pt.data$subj==i])) ) }
N <- length(bx.subj)

bx.full <- as.data.frame(cbind(bx.subj, bx.time))
names(bx.full) <- c("subj", "time.int")


#add numeric date and age at start of interval and at time of biopsy; add whether bx was performed during interval; add bx findings
bx.full$int.dt.num <- bx.full$int.age <- vector(length=N)
bx.full$bx.time <- bx.full$bx.dt.num <- bx.full$bx.age <- vector(length=N)
bx.full$bx.here <- bx.full$num.bx<- vector(length=N) ## bx.here indicates whether there was biopsy occured in that time interval 111422 Zitong
bx.full$rc <- bx.full$pgg <-  vector(length=N)
bx.full$bx.time.min <- bx.full$bx.dt.num.min <- bx.full$bx.age.min <- rep(NA, N) #for intervals with 2 bx
bx.full$npc <- bx.full$mpc <- bx.full$ncs <- bx.full$mri <- vector(length=N)
bx.full$npc.min <- bx.full$mpc.min <- bx.full$ncs.min <- bx.full$mri.min <- vector(length=N)
bx.full$lat <- bx.full$lat.min <- vector(length=N)
bx.full$pirads <- vector(length=N) #1126 Zitong

#data at time=0
bx.full$bx.time[bx.full$time.int==0] <- 0
bx.full$bx.here[bx.full$time.int==0] <- 1
bx.full$rc[bx.full$time.int==0] <- 0
bx.full$pgg[bx.full$time.int==0] <- 1

for(i in 1:n){
	bx.full$npc[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$npc[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"]
	bx.full$mpc[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$mpc[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"]
	bx.full$ncs[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$ncs[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"]
	bx.full$mri[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$mri[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"]
	bx.full$lat[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$lat[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"]
  if(date.pull == as.numeric(as.Date("2022-11-07"))){ #if new data
    bx.full$pirads[bx.full$subj==pt.data$subj[i] & bx.full$time.int==0] <- bx.data$max_pirads[bx.data$subj==pt.data$subj[i] & bx.data$Which_Biopsy=="D"] #1126 Zitong
  }
}

#variable values that do not depend on patient biopsy data, just pt.data
for(i in 1:n){
#at time=0
	bx.full$bx.dt.num[bx.full$time.int==0 & bx.full$subj==pt.data$subj[i]] <- pt.data$dx.dt.num[i]
	bx.full$bx.age[bx.full$time.int==0 & bx.full$subj==pt.data$subj[i]] <- pt.data$dx.age[i]
#at all times
	bx.full$int.dt.num[bx.full$subj==pt.data$subj[i]] <- pt.data$dx.dt.num[i] + bx.full$time.int[bx.full$subj==pt.data$subj[i]]*365
	bx.full$int.age[bx.full$subj==pt.data$subj[i]] <- pt.data$dx.age[i] + bx.full$time.int[bx.full$subj==pt.data$subj[i]]}


#variable values that do depend on patient biopsy data
options(warn = 0)
for(j in 1:N){
	if(bx.full$time.int[j]>0){ #post-dx biopsies
		bx.data$use<-rep(0,n_bx) #clearing existing values in this variable
		bx.data$use[ bx.data$subj==bx.full$subj[j] & bx.data$time.since.dx>(bx.full$time.int[j]-1) & bx.data$time.since.dx<=bx.full$time.int[j] ] <- 1 #identifying biopsies for this patient in the time interval of interest


		if(sum(bx.data$use) > 0){ #if this patient had any biopsies in this time interval
		  bx.full$bx.time[j] <- max(bx.data$time.since.dx[bx.data$use==1]) #want to use max in case patient had multiple biopies in an interval
			bx.full$bx.dt.num[j] <- max(bx.data$bx.dt.num[bx.data$use==1])
			bx.full$bx.age[j] <- bx.full$bx.time[j] + pt.data$dx.age[pt.data$subj==bx.full$subj[j]]
			bx.full$bx.here[j] <- 1
			bx.full$rc[j] <- bx.data$rc[bx.data$use==1 & bx.data$bx.dt.num==bx.full$bx.dt.num[j]]
			bx.full$pgg[j] <- bx.data$pgg[bx.data$use==1 & bx.data$bx.dt.num==bx.full$bx.dt.num[j]]
			bx.full$num.bx[j] <- sum(bx.data$use)
      #if(sum(bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j])>1){print(j)}
						bx.full$npc[j] <- bx.data$npc[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]]
						bx.full$mpc[j] <- bx.data$mpc[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]]
						bx.full$ncs[j] <- bx.data$ncs[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]]
						bx.full$mri[j] <- bx.data$mri[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]]
						bx.full$lat[j] <- bx.data$lat[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]]
						if(date.pull == as.numeric(as.Date("2022-11-07"))){ #if new data
						bx.full$pirads[j] <- bx.data$max_pirads[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time[j]] ##1126 Zitong
						}
			if(bx.full$num.bx[j]>1){ #patients with multiple biopsies
				bx.full$bx.time.min[j] <- min(bx.data$time.since.dx[bx.data$use==1])
				bx.full$bx.dt.num.min[j] <- min(bx.data$bx.dt.num[bx.data$use==1])
				bx.full$bx.age.min[j] <- bx.full$bx.time.min[j] + pt.data$dx.age[pt.data$subj==bx.full$subj[j]]

				bx.full$npc.min[j] <- bx.data$npc[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]
				bx.full$mpc.min[j] <- bx.data$mpc[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]
				bx.full$ncs.min[j] <- bx.data$ncs[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]
				bx.full$mri.min[j] <- bx.data$mri[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]
				bx.full$lat.min[j] <- bx.data$lat[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]
				if(date.pull == as.numeric(as.Date("2022-11-07"))){ #if new data
				bx.full$pirads.min[j] <- bx.data$max_pirads[bx.data$use==1 & bx.data$time.since.dx==bx.full$bx.time.min[j]]##1126 Zitong
				}
			} #we know that rc=0 at this bx bc all bx after rc were removed from dataset

		}

		else{ #if the patient didn't get any biopsies in this interval
			bx.full$bx.time[j] <- bx.full$bx.dt.num[j] <- bx.full$bx.age[j] <-
			  bx.full$rc[j] <- bx.full$pgg[j] <-
			  bx.full$npc[j] <- bx.full$mpc[j] <- bx.full$npc.min[j] <- bx.full$mpc.min[j] <-
			  bx.full$ncs[j] <- bx.full$mri[j] <- bx.full$ncs.min[j] <- bx.full$mri.min[j] <-
			  bx.full$lat[j] <- bx.full$lat.min[j] <- NA
			  if(date.pull == as.numeric(as.Date("2022-11-07"))){ #if new data
			    bx.full$pirads[j] <- bx.full$pirads.min[j] <-  NA ##1126 Zitong
		  	}
			bx.full$bx.here[j] <- bx.full$num.bx[j] <- 0 }	}
} #

#how many biopsies in this interval?
if(max(bx.full$num.bx)>2){
	warning(paste0(length(unique(bx.full$subj[bx.full$num.bx>2])), " patients have at least one annual interval with more than 2 biopsies. Only data from two biopsies (the first and last in the annual interval) will be used."))}

#this checks whether code for laterality assignments captured all bilateral findings
if(!sum(bx.data$lat, na.rm=T)==(sum(bx.full$lat, na.rm=T) + sum(bx.full$lat.min, na.rm=T))){
  warning("Number of biopsies with bilateral PCa first found not equal.")}


#did patient get survery in this interval?
bx.full$surg<-rep(0,N)
for(i in 1:n){
	if(pt.data$surg[i]==1){
	  bx.full$surg[bx.full$subj==pt.data$subj[i] & bx.full$time.int==ceiling(pt.data$time.until.tx[i])] <- 1}}

#did patient have treatment in this interval?
bx.full$tx<-rep(0,N)
for(i in 1:n){
  if(!is.na(pt.data$time.until.tx[i])){
    bx.full$tx[bx.full$subj==pt.data$subj[i] & bx.full$time.int==ceiling(pt.data$time.until.tx[i])] <- 1}}

## add prostate volume
bx.full$avg.vol<-bx.full$std.vol <- vector(length=N)
for(i in 1:n){
  bx.full$std.vol[bx.full$subj==pt.data$subj[i]] <- pt.data$std.vol[i]
  bx.full$avg.vol[bx.full$subj==pt.data$subj[i]] <- pt.data$avg.vol[i] }


#predictors of biopsy and surgery received (looking back at all biopsy findings that would influence decisions at beginning and ending of annual interval)
#number of previous biopsies
bx.full$num.prev.bx.start <- bx.full$num.prev.bx.end <- vector(length=N)

bx.full$num.prev.bx.start[bx.full$time.int==0] <- 0
bx.full$num.prev.bx.end[bx.full$time.int==0] <- 1

#did patient experience reclssification earlier (including this interval)?
bx.full$prev.rc <- vector(length=N)
#RC with biopsy grade >=3 or >=4
bx.full$prev.grade.rc.3 <- bx.full$prev.grade.rc.4 <- vector(length=N)


#max previous npc and mpc
bx.full$max.prev.npc.start <- bx.full$max.prev.npc.end <- vector(length=N)
bx.full$max.prev.mpc.start <- bx.full$max.prev.mpc.end <- vector(length=N)
bx.full$cuml.prev.ncs.start <- bx.full$cuml.prev.ncs.end <- vector(length=N)
bx.full$cuml.prev.mri.start <- bx.full$cuml.prev.mri.end <- vector(length=N)
bx.full$prev.lat.start <- bx.full$prev.lat.end <- vector(length=N)


#for loop to populate # prev bx and prev results fields
for(j in 1:N){if(bx.full$time.int[j]>0){

	bx.full$num.prev.bx.start[j] <-
	  sum(bx.full$num.bx[bx.full$subj==bx.full$subj[j]
	                     & bx.full$time.int<bx.full$time.int[j]])
	bx.full$num.prev.bx.end[j] <-
	  sum(bx.full$num.bx[bx.full$subj==bx.full$subj[j]
	                     & bx.full$time.int<=bx.full$time.int[j]])

	bx.full$prev.rc[j] <- max(bx.full$rc[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]], na.rm=T)

	if(bx.full$prev.rc[j]==1){
		bx.full$prev.grade.rc.3[j] <- as.numeric(max(bx.full$pgg[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]], na.rm=T)>=3)
		bx.full$prev.grade.rc.4[j] <- as.numeric(max(bx.full$pgg[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]], na.rm=T)==4) }
	if(bx.full$prev.rc[j]==0){
		bx.full$prev.grade.rc.3[j] <- bx.full$prev.grade.rc.4[j] <- 0 }

	bx.full$max.prev.npc.start [j] <-
	  max(c(bx.full$npc[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]],
	        bx.full$npc.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]]), na.rm=T)
	bx.full$max.prev.npc.end [j] <-
	  max(c(bx.full$npc[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]],
	        bx.full$npc.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]]), na.rm=T)

	bx.full$max.prev.mpc.start [j] <-
	  max(c(bx.full$mpc[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]],
	        bx.full$mpc.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]]), na.rm=T)
	bx.full$max.prev.mpc.end [j] <-
	  max(c(bx.full$mpc[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]],
	        bx.full$mpc.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]]), na.rm=T)

	bx.full$cuml.prev.ncs.start [j] <-
	  sum(c(bx.full$ncs[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]],
	        bx.full$ncs.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]]), na.rm=T)
	bx.full$cuml.prev.ncs.end [j] <-
	  sum(c(bx.full$ncs[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]],
	        bx.full$ncs.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]]), na.rm=T)

	bx.full$cuml.prev.mri.start [j] <-
	  sum(c(bx.full$mri[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]],
	        bx.full$mri.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]]), na.rm=T)
	bx.full$cuml.prev.mri.end [j] <-
	  sum(c(bx.full$mri[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]],
	        bx.full$mri.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]]), na.rm=T)

	bx.full$prev.lat.start [j] <-
	  max(c(bx.full$lat[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]],
	        bx.full$lat.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<bx.full$time.int[j]]), na.rm=T)
	bx.full$prev.lat.end [j] <-
	  max(c(bx.full$lat[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]],
	        bx.full$lat.min[bx.full$subj==bx.full$subj[j] & bx.full$bx.here==1 & bx.full$time.int<=bx.full$time.int[j]]), na.rm=T)

} }


bx.full$npc.min[is.na(bx.full$bx.dt.num.min)]<-NA
bx.full$mpc.min[is.na(bx.full$bx.dt.num.min)]<-NA
bx.full$lat.min[is.na(bx.full$bx.dt.num.min)]<-NA

bx.full$freq.bx.start <- bx.full$num.prev.bx.start /bx.full$time.int
bx.full$freq.bx.end <- bx.full$num.prev.bx.end /(bx.full$time.int+1)

bx.full$num.npc0.start<-bx.full$num.npc0.end<-vector(length=N)
bx.full$prop.npc0.start<-bx.full$prop.npc0.end<-vector(length=N)


for(j in 1:N){
  if(bx.full$time.int[j]>0){
    subj_j<-bx.full$subj[j]
    time_j<-bx.full$time.int[j]

    bx.full$num.npc0.start[j] <-
      sum(bx.full$npc==0 & !is.na(bx.full$npc)
          & bx.full$bx.here==1 & !is.na(bx.full$bx.here)
          & bx.full$subj==subj_j & bx.full$time.int<time_j) +
      sum(bx.full$npc.min==0 & !is.na(bx.full$npc.min)
          & bx.full$num.bx > 1
          & bx.full$subj==subj_j & bx.full$time.int<time_j)

    bx.full$num.npc0.end[j] <-
      sum(bx.full$npc==0 & !is.na(bx.full$npc)
          & bx.full$bx.here==1 & !is.na(bx.full$bx.here)
          & bx.full$subj==subj_j & bx.full$time.int<=time_j) +
      sum(bx.full$npc.min==0 & !is.na(bx.full$npc.min)
           & bx.full$num.bx > 1
          & bx.full$subj==subj_j & bx.full$time.int<=time_j)

    bx.full$prop.npc0.start[j] <-
      bx.full$num.npc0.start[j]/bx.full$num.prev.bx.start[j]

    bx.full$prop.npc0.end[j] <-
      bx.full$num.npc0.end[j]/bx.full$num.prev.bx.end[j]
      }
}



#did volume reclassification occur? may also influence decisions
bx.full$prev.vol.rc.start <- as.numeric(bx.full$max.prev.npc.start>2 | bx.full$max.prev.mpc.start>50)
bx.full$prev.vol.rc.end <- as.numeric(bx.full$max.prev.npc.end>2 | bx.full$max.prev.mpc.end>50)




#Indicate that patients aren't eligible for bx post-RC.
#Actual biopsy information on post-RC biopsies already removed from bx.data
for(i in 1:n){
	if(pt.data$rc[i]==1){
		rc.time<-bx.full$time.int[bx.full$subj==pt.data$subj[i] & bx.full$rc==1 & !is.na(bx.full$rc)]
		bx.full$bx.here[bx.full$subj==pt.data$subj[i] & bx.full$time.int>rc.time]<-NA}}

data.check(condition=as.logical(length(unique(bx.full$subj[is.na(bx.full$bx.here)])) == (sum(pt.data$rc==1 & is.na(pt.data$time.until.tx) & !ceiling(pt.data$time.until.rc)==ceiling(pt.data$time.fup)) + sum(!is.na(pt.data$time.until.rc) & !is.na(pt.data$time.until.tx) & !ceiling(pt.data$time.until.rc)==ceiling(pt.data$time.until.tx) ) ) ), message="Error in defining biopsy intervals. Email Yates; she will check orginal script.")



#present PSA at the beginning and end of interval (for biopsy and surgery decisions, respectively)
bx.full$psa.pres.bx <- bx.full$psa.pres.surg <- rep(NA,N)
for(j in 1:N){
	if(bx.full$time.int[j]>0){
	psa.use.bx <- psa.use.surg <- date.use <- NA
	if(bx.full$bx.here[j]==1 & !is.na(bx.full$bx.here[j])){
#patients with biopsies in the interval
		if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] & psa.data$time.since.dx>bx.full$time.int[j]-0.5) >0 ){ #with at least one PSA in prior window
#for biopsy
			psa.use.bx <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] & psa.data$time.since.dx>bx.full$time.int[j]-0.5]
			bx.full$psa.pres.bx[j] <- mean(psa.use.bx)}
		if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] + (11/12) & psa.data$time.since.dx>bx.full$time.int[j]-0.5) >0 ){
#for surgery
			psa.use.surg <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j]+11/12 & psa.data$time.since.dx>bx.full$time.int[j]-0.5]
			bx.full$psa.pres.surg[j] <- mean(psa.use.surg) }
			}

#for patients without biopsies in the interval
		else{
			if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-0.5) > 0){
			psa.use.bx <- psa.use.surg <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-0.5]
			bx.full$psa.pres.bx[j] <- bx.full$psa.pres.surg[j] <- mean(psa.use.bx)} }

		if(sum(is.na(psa.use.bx))>0){
			if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j])>0){
			date.use <- max(psa.data$dt.num[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j]])
			if(sum(is.na(psa.use.surg))>0){
				bx.full$psa.pres.bx[j] <- bx.full$psa.pres.surg[j] <- unique(psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$dt.num==date.use])} ## add unique(), RHS returns 2 values
				else{
					bx.full$psa.pres.bx[j]  <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$dt.num==date.use]
				}
			} }
			} }

summary(bx.full$psa.pres.bx[!bx.full$time.int==0])
summary(bx.full$psa.pres.surg[!bx.full$time.int==0])
#hist(bx.full$psa.pres.surg[!bx.full$time.int==0])

bx.full$psa.pres.bx[is.na(bx.full$psa.pres.bx)]<-mean(bx.full$psa.pres.bx,na.rm=T)
bx.full$psa.pres.surg[is.na(bx.full$psa.pres.surg)]<-mean(bx.full$psa.pres.surg,na.rm=T)


#change in PSA ("traj"=trajectory)
bx.full$psa.traj.bx <- bx.full$psa.traj.surg <- rep(NA,N)

for(j in 1:N){
	if(bx.full$time.int[j]>0){
	psa.use.bx <- psa.use.surg <- sum.use <- times.use <- NA

#patients with biopsies
	if(bx.full$bx.here[j]==1 & !is.na(bx.full$bx.here[j])){

		if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] & psa.data$time.since.dx>bx.full$time.int[j]-2) >1 ){

			psa.use.bx <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] & psa.data$time.since.dx>bx.full$time.int[j]-2]
			times.use <- psa.data$time.since.dx[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$bx.time[j] & psa.data$time.since.dx>bx.full$time.int[j]-2]
			bx.full$psa.traj.bx[j] <- lm(psa.use.bx~times.use)$coefficients[2] }

		if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j]+(11/12) & psa.data$time.since.dx>bx.full$time.int[j]-2) > 1 ){

			psa.use.surg <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-2]
			times.use <- psa.data$time.since.dx[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-2]
			bx.full$psa.traj.surg[j] <- lm(psa.use.surg~times.use)$coefficients[2] } }


#patients without biopsies in interval
		else{
			if(sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-2) > 1){
			psa.use.bx <- psa.use.surg <- psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-2]
			times.use <- psa.data$time.since.dx[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<(bx.full$time.int[j]+(11/12)) & psa.data$time.since.dx>bx.full$time.int[j]-2]
			bx.full$psa.traj.bx[j] <- bx.full$psa.traj.surg[j] <- lm(psa.use.bx~times.use)$coefficients[2] } }

		if(sum(is.na(psa.use.bx))>0){
			sum.use <- sum(psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j])
			if(sum.use>1){
			times.use <- sort(psa.data$time.since.dx[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx<bx.full$time.int[j]])[(sum.use-1):sum.use]
			psa.use.bx <- c(mean(psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx==times.use[1]]), mean(psa.data$log.psa[psa.data$subj==bx.full$subj[j] & psa.data$time.since.dx==times.use[2]]) )
				if(sum(is.na(psa.use.surg))>0){
					bx.full$psa.traj.bx[j] <- bx.full$psa.traj.surg[j] <- lm(psa.use.bx~times.use)$coefficients[2] }
					else{ bx.full$psa.traj.bx[j] <- lm(psa.use.bx~times.use)$coefficients[2] } 	}
} }
			}

#summary(bx.full$psa.traj.surg[!bx.full$time.int==0])
#hist(bx.full$psa.traj.surg[!bx.full$time.int==0], breaks=20)

#quantile(bx.full$psa.traj.bx, p=c(0.025, 0.975), na.rm=T)
bx.full$psa.traj.bx[bx.full$psa.traj.bx <
                      quantile(bx.full$psa.traj.bx, p=0.025, na.rm=T)] <-
  quantile(bx.full$psa.traj.bx, p=0.025, na.rm=T)

bx.full$psa.traj.surg[bx.full$psa.traj.bx <
                      quantile(bx.full$psa.traj.surg, p=0.025, na.rm=T)] <-
  quantile(bx.full$psa.traj.surg, p=0.025, na.rm=T)


bx.full$psa.traj.bx[bx.full$psa.traj.bx >
                      quantile(bx.full$psa.traj.bx, p=0.975, na.rm=T)] <-
  quantile(bx.full$psa.traj.bx, p=0.975, na.rm=T)

bx.full$psa.traj.surg[bx.full$psa.traj.bx >
                        quantile(bx.full$psa.traj.surg, p=0.975, na.rm=T)] <-
  quantile(bx.full$psa.traj.surg, p=0.975, na.rm=T)


bx.full$psa.traj.bx[is.na(bx.full$psa.traj.bx)]<-0
bx.full$psa.traj.surg[is.na(bx.full$psa.traj.surg)]<-0


### 10. add MRI data with MRI findings --------------
###Zitong 113022
mri.data <- read_csv(paste(location.of.data, name.of.mri.data, sep = "/"))

if(date.pull == as.numeric(as.Date("2022-11-07"))){
mri.find.data <- read_csv(paste(location.of.data, name.of.mri.findings.data, sep = "/"))
(mri.prob <- problems(x=.Last.value))
mri.data[mri.prob$row-1, mri.prob$col] <- as.numeric(gsub("\\.00|,", "", mri.prob$actual))
mri.data.new <- mri.data %>% left_join(mri.find.data, by = "mri_id")
#assign all missing pirads but no_measurable_disease=1 to have pirads 1
mri.data.new<- mri.data.new %>% mutate(pirads_update  = 
                                         ifelse(is.na(pirads) & no_measurable_disease == 1, 
                                                1,
                                                pirads))
}else{
  mri.data.new <- mri.data
}


# ### 11. Save data
save(pt.data, psa.data, bx.full,mri.data.new,
     file=paste0(location.of.generated.files,"/IOP-data-shaping-work-space-6.15-withMRI.RData"))

#clean-up workspace
rm(bx.data, bx.subj, bx.time, d,  demo.data, i, j, max.fup, max.time.until.surg, n, N, n_psa, n_bx, n_tx, ordered, psa.use.bx, psa.use.surg, rc.time, rm.id, sum.use, times.use)

