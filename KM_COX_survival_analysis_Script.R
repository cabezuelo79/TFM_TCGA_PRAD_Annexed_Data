##############################################################
#######        KAPLAN-MEIER SURVIVAL ANALYSIS        #########        
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(survival)
library(gtsummary)

# read PRAD_RNA file 
rna <- read.table('RNA/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t') 
#samples not paired must be removed in KM paired analysis

# take off first row cause don't need it
rna <- rna[-1,]

# read the Clinical file, transposed it to keep the clinical feature title as column name
clinical <- t(read.table('Clinical/PRAD.merged_only_clinical_clin_format.txt',header=T, row.names=1, sep='\t'))
#samples not paired must be removed in KM paired analysis

# remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
 x <- as.matrix(x)
 x <- t(apply(x,1,as.numeric))
 r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
 remove <- which(r > dim(x)[2]*0.5)
 return(remove)
}
 remove <- rem(rna)
 rna <- rna[-remove,]

# see the values
table(substr(colnames(rna),14,14))

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')

# apply voom function from limma package to normalize the data
library('limma')
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% t_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}

rna_vm  <- vm(rna)
colnames(rna_vm) <- gsub('\\.','-',substr(colnames(rna),1,12))

# check data 
hist(rna_vm)

# remove the old "rna" 
rm(rna)

### z = [(value gene X in tumor Y)-(mean gene X in normal)]/(standard deviation X in normal) ###

# calculate z-scores
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
z_rna <- scal(rna_vm[,t_index],rna_vm[,n_index])


# set the rownames keeping only gene name
rownames(z_rna) <- sapply(rownames(z_rna), function(x) unlist(strsplit(x,'\\|'))[[1]])

rm(rna_vm) #we don't need it anymore

# match the patient ID in clinical data with the colnames of z_rna
clinical <- as.data.frame(clinical)
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(clinical$IDs %in% colnames(z_rna)) 

# get the columns that contain data we can use
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))

# collapse follow ups together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if (sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and gleason
ind_keep <- grep('patient.stage_event.gleason_grading.gleason_score',colnames(clinical))
gleason <- as.matrix(clinical[,ind_keep])
gleason_collapsed <- c()
for (i in 1:dim(gleason)[1]){
  if (sum (is.na(gleason[i,])) < dim(gleason)[2]){
    m <- max(gleason[i,],na.rm=T)
    gleason_collapsed <- c(gleason_collapsed,m)
  } else {
    gleason_collapsed <- c(gleason_collapsed,'NA')
  }
}

# and PSA
ind_keep <- grep('patient.stage_event.psa.psa_value',colnames(clinical))
psa <- as.matrix(clinical[,ind_keep])
psa_collapsed <- c()
for (i in 1:dim(psa)[1]){
  if (sum (is.na(psa[i,])) < dim(psa)[2]){
    m <- max(psa[i,],na.rm=T)
    psa_collapsed <- c(psa_collapsed,m)
  } else {
    psa_collapsed <- c(psa_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed,gleason_collapsed,psa_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days','gleason_score','PSA')


# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}

# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}

# create a new vector time by merging the two previous ones
all_clin$new_time_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_time)))){
  all_clin$new_time_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_time))[i]),
                                   as.numeric(as.character(all_clin$new_death))[i],as.numeric(as.character(all_clin$new_time))[i])
}

# create vector for tumour status censoring
table(clinical$patient.person_neoplasm_cancer_status)

all_clin$tumour_event <- ifelse(clinical$patient.person_neoplasm_cancer_status == 'tumor free', 0,1)

# create vector for aggressivity censoring based on gleason 
table(clinical$patient.stage_event.gleason_grading.gleason_score)

all_clin$aggressivity <- ifelse(clinical$patient.stage_event.gleason_grading.gleason_score == '6', 0,1) # The lowest gleason value in samples is 6, above that is considered aggressive

# create vector for PSA censoring
table(clinical$patient.stage_event.psa.psa_value) # the values are input without an uniform format so it is decided to introduce them manually. DELETE NON PAIRED SAMPLES WHEN PERFORM PAIRED ANALYSIS!!!

# values under 10 ng/mL are considered '0' values over 10 are considered '1' 
all_clin$psa_score <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,NA,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,NA,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,NA,NA,NA,0,1,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,NA,0,0,0,0,0,NA,0,0,1,0,0,0,NA,NA,NA,0,0,NA,0,0,0,NA,0,0,0,0,0,0,0,0,NA,0,0,NA,0,NA,0,0,0,0,NA,0,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,1,0,0,0,0,NA,NA,NA,NA,NA,0,NA,NA,0,0,0,0,0,0,0,NA,1,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,NA,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,NA,NA,NA,0,NA,NA,0,0,0,NA,0,NA,0,0,0,0,0,0,0,0,1,0,0,0)

#finally add row.names to clinical
rownames(all_clin) <- clinical$IDs

# create event vector for RNASeq data
event_rna <- t(apply(z_rna, 1, function(x) ifelse(abs(x) > 1.96,1,0)))

# since we need the same number of patients in both clinical and RNASeq data take the indices for the matching samples
ind_tum <- which(unique(colnames(z_rna)) %in% rownames(all_clin))
ind_clin <- which(rownames(all_clin) %in% colnames(z_rna))

# pick gene of interest
ind_gene <- which(rownames(z_rna) == 'KCNG1')

# check how many altered samples we have
table(event_rna[ind_gene,])



############## run survival analysis for TUMOR_STATUS ##############

s <- survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$tumour_event[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$tumour_event[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

pv

# plot the data
plot(survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$tumour_event[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('TUMOUR STATUS: PRESENT',rownames(z_rna)[ind_gene],sep='\n'))

# add lines for the median 
x1 <- ifelse ( is.na(as.numeric(summary(s)$table[,'median'][1])),'NA',as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if (x1 != 'NA' & x2 != 'NA'){
  lines(c(0,x1),c(0.5,0.5),col='blue')
  lines(c(x1,x1),c(0,0.5),col='black')
  lines(c(x2,x2),c(0,0.5),col='red')
}

# add legend
legend(x='topright', 1800,0.995,legend=paste('p.value = ',pv[[1]],sep=''),bty='n',cex=1.4)
legend(x='right', max(as.numeric(as.character(all_clin$aggressivity)[ind_clin]),na.rm = T)*0.7,0.94,
       legend=c(paste('Tumour absent=',x1),paste('Tumour present',x2)),bty='n',cex=1.3,lwd=3,col=c('black','red'))


# COX HAZARD ESTIMATION
coxph(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$tumour_event[ind_clin])~event_rna[ind_gene,ind_tum])%>%gtsummary::tbl_regression(exp = TRUE) 



############## run survival analysis for TUMOR AGGRESSIVITY ###############

s <- survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$aggressivity[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$aggressivity[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

pv

# plot the data
plot(survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$aggressivity[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('aggressivity',rownames(z_rna)[ind_gene],sep='\n'))

# add lines for the median 
x1 <- ifelse ( is.na(as.numeric(summary(s)$table[,'median'][1])),'NA',as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if (x1 != 'NA' & x2 != 'NA'){
  lines(c(0,x1),c(0.5,0.5),col='blue')
  lines(c(x1,x1),c(0,0.5),col='black')
  lines(c(x2,x2),c(0,0.5),col='red')
}

# add legend
legend(x='topright', 1800,0.995,legend=paste('p.value = ',pv[[1]],sep=''),bty='n',cex=1.4)
legend(x='right', max(as.numeric(as.character(all_clin$aggressivity)[ind_clin]),na.rm = T)*0.7,0.94,
       legend=c(paste('Low Gleason score=',x1),paste('High Gleason score=',x2)),bty='n',cex=1.3,lwd=3,col=c('black','red'))


#COX HAZARD ESTIMATION
coxph(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$aggressivity[ind_clin])~event_rna[ind_gene,ind_tum])%>%gtsummary::tbl_regression(exp = TRUE) 


##############    run survival analysis for PSA LEVELS    ###############

s <- survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum])
s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))

# extraect the p.value
pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]

pv

# plot the data
plot(survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum]),
     col=c(1:3), frame=F, lwd=2,main=paste('PSA_LEVELS',rownames(z_rna)[ind_gene],sep='\n'))

# add lines for the median 
x1 <- ifelse ( is.na(as.numeric(summary(s)$table[,'median'][1])),'NA',as.numeric(summary(s)$table[,'median'][1]))
x2 <- as.numeric(summary(s)$table[,'median'][2])
if (x1 != 'NA' & x2 != 'NA'){
  lines(c(0,x1),c(0.5,0.5),col='blue')
  lines(c(x1,x1),c(0,0.5),col='black')
  lines(c(x2,x2),c(0,0.5),col='red')
}

# add legend
legend(x='topright', 1800,0.995,legend=paste('p.value = ',pv[[1]],sep=''),bty='n',cex=1.4)
legend(x='right', max(as.numeric(as.character(all_clin$psa_score)[ind_clin]),na.rm = T)*0.7,0.94,
       legend=c(paste('NotAltered=',x1),paste('Altered=',x2)),bty='n',cex=1.3,lwd=3,col=c('black','red'))


# COX HAZARD ESTIMATION
coxph(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum])%>%gtsummary::tbl_regression(exp = TRUE) 


##############################################################
#######               TESTING A GENE SET             #########        
##############################################################

# ALL THE PREVIOUS STEPS UNTIL 'pick gene of interest' ARE COMMON

gene_set <- #VECTOR WITH SET OF DEG SYMBOLS

counter = 0

all_pv <- data.frame()
all_gene <- data.frame()

for (gene in (gene_set)){
  counter=counter +1

  # pick your gene of interest
  ind_gene <- which(rownames(z_rna) == gene)
  
  # check how many altered samples we have
  table(event_rna[ind_gene,])
  
  # run survival analysis for PSA SCORE as example 
  s <- survfit(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum])
  s1 <- tryCatch(survdiff(Surv(as.numeric(as.character(all_clin$new_time_death))[ind_clin],all_clin$psa_score[ind_clin])~event_rna[ind_gene,ind_tum]), error = function(e) return(NA))
  
  # extract the p.value
  pv <- ifelse ( is.na(s1),next,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),3)))[[1]]
  all_pv[counter,'pv']<-pv
  all_gene[counter,'gene']<-gene  
  
}

gene_pv <- data.frame(all_gene,all_pv)

gene_pv
