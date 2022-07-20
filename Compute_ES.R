library("tidyverse")
library("dplyr")
library("readxl")
library("writexl")

directory = "/Users/shreygrover/Desktop/MetaAnalysis_SciTranslMed/Sci_Transl_Med_Submission_2/Analysis/Scripts_Revised/Meta_Datafiles_2021"

# functions

# computes Hedges' G for within-subjects experiments
compute_G_within <- function(n, M1, M2, SD1, SD2, r)
{
  J = 1 - 3/(4*(n-1)-1) # correction factor J 
  d = ((M1 - M2) / (sqrt(SD1^2 + SD2^2 - 2*r*SD1*SD2)))*sqrt(2*(1-r)) # Cohen's drm
  Vd = (1/n + d^2/(2*n))*2*(1-r) # variance of drm  
  g = J*d # G
  Vg = J^2*Vd # variance of G
  SEg = sqrt(Vg) # SE of G
  return(list(g,SEg))
}

# computes Hedges' G for between-subjects experiments
compute_G_between <- function(n1, n2, M1, M2, SD1, SD2)
{
  J = 1 - 3/(4*(n1+n2-2)-1) # correction factor J 
  d = (M1 - M2) / (sqrt(((n1-1)*SD1^2 + (n2-1)*SD2^2)/(n1+n2-2))) # Cohen's d
  Vd = ((n1+n2)/(n1*n2)) + ((d^2)/(2*(n1+n2))) # variance of d  
  g = J*d # G
  Vg = (J^2)*Vd # variance of G
  SEg = sqrt(Vg) # SE of G
  return(list(g,SEg))
}

# computes Hedges' G for within-subjects experiments using difference score and its SD
compute_G_diffScore <- function(n, Mdiff, SDdiff, r)
{
  J = 1 - 3/(4*(n-1)-1) # correction factor J 
  d = (Mdiff/SDdiff)*sqrt(2*(1-r)) # Cohen's drm
  Vd = (1/n + d^2/(2*n))*2*(1-r) # variance of drm 
  g = J*d # G
  Vg = (J^2)*Vd # variance of G
  SEg = sqrt(Vg) # SE of G
  return(list(g,SEg))
}

# computes Hedges' G for between-subjects experiments using independent samples t-test statistic
compute_G_indT <- function(n1,n2,t)
{
  J = 1 - 3/(4*(n1+n2-2)-1) # correction factor J 
  d = t*sqrt(1/n1+1/n2) # Cohen's d
  Vd = ((n1+n2)/(n1*n2)) + ((d^2)/(2*(n1+n2))) # variance of d  
  g = J*d # G
  Vg = (J^2)*Vd # variance of G
  SEg = sqrt(Vg) # SE of G
  return(list(g,SEg))
}

# Computation of Hedges' G

# choose a dataset (outcome, hypothesis)

# metadata <- read_excel(paste0(directory, "/ES_rawdata/ES_raw_hypothesis.xlsx"))
# metadata <- read_excel(paste0(directory, "/ES_rawdata/ES_raw_outcome.xlsx"))
metadata <- read_excel(paste0(directory, "/ES_rawdata/ES_raw_NEW_EFFECTS_nov21.xlsx"))

# convert G and seG columns to numberic
metadata$G <- as.numeric(metadata$G)
metadata$seG <- as.numeric(metadata$seG)
attach(metadata)

# correlation between condition in within-subjects experiments
r = 0.5 # default 0.5; sensitivity: 0.3, 0.7
metadata$corr_coef[is.na(metadata$corr_coef)] <- r # assign only to studies where we don't know r

# Compute G using depending on study design
  metadata[metadata$ES_method=='WITHIN',c('G','seG')]<-compute_G_within(
  n   = n[metadata$ES_method=='WITHIN'],
  M1  = M1[metadata$ES_method=='WITHIN'],
  M2  = M2[metadata$ES_method=='WITHIN'],
  SD1 = SD1[metadata$ES_method=='WITHIN'],
  SD2 = SD2[metadata$ES_method=='WITHIN'],
  r=r
)

metadata[metadata$ES_method=='BETWEEN',c('G','seG')]<-compute_G_between(
  n1  = n1[metadata$ES_method=='BETWEEN'],
  n2  = n2[metadata$ES_method=='BETWEEN'],
  M1  = M1[metadata$ES_method=='BETWEEN'],
  M2  = M2[metadata$ES_method=='BETWEEN'],
  SD1 = SD1[metadata$ES_method=='BETWEEN'],
  SD2 = SD2[metadata$ES_method=='BETWEEN']
)

metadata[metadata$ES_method=='DIFF_SCORE',c('G','seG')]<-compute_G_diffScore(
  n      = n[metadata$ES_method=='DIFF_SCORE'],
  Mdiff  = M1[metadata$ES_method=='DIFF_SCORE'],
  SDdiff = SD1[metadata$ES_method=='DIFF_SCORE'],
  r      = r
)

metadata[metadata$ES_method=='IND_t_test',c('G','seG')]<-compute_G_indT(
  n1 = n1[metadata$ES_method=='IND_t_test'],
  n2 = n2[metadata$ES_method=='IND_t_test'],
  t  = t[metadata$ES_method=='IND_t_test']
)

metadata$meanG <- metadata$G
metadata$meanVarG <- metadata$seG^2

# Combine effects

# array of markers for combining effects
es_markers <- sort(unique(metadata$comb_ef))
es_markers[1] <- NaN 

for(i in 1:(length(es_markers)-1)){
  # take the mean of G of effects that are combined
  metadata[metadata$comb_ef == es_markers[i+1],]$meanG <- mean(subset(metadata, comb_ef==es_markers[i+1])$G)
  temp <- subset(metadata, comb_ef==es_markers[i+1]) # temporary subset of data
  n_exp <- length(temp$experiment) # number of experiments to combine
  r = 0.5 # correlation between effects
  
  # compute mean SDs of combined effects
  if(n_exp==5) {
    metadata[metadata$comb_ef == es_markers[i+1],]$meanVarG <- (1/n_exp)^2 * (sum(temp$seG^2, na.rm = TRUE) + 
                                 2*r*temp$seG[1]*temp$seG[2] + 2*r*temp$seG[1]*temp$seG[3] + 2*r*temp$seG[1]*temp$seG[4] + 
                                 2*r*temp$seG[1]*temp$seG[5] + 2*r*temp$seG[2]*temp$seG[3] + 2*r*temp$seG[2]*temp$seG[4] + 
                                 2*r*temp$seG[2]*temp$seG[5] + 2*r*temp$seG[3]*temp$seG[4] + 2*r*temp$seG[3]*temp$seG[5] +
                                 2*r*temp$seG[4]*temp$seG[5])   
  } else if(n_exp==4) {
    metadata[metadata$comb_ef == es_markers[i+1],]$meanVarG <- (1/n_exp)^2 * (sum(temp$seG^2, na.rm = TRUE) + 
                                 2*r*temp$seG[1]*temp$seG[2] + 2*r*temp$seG[1]*temp$seG[3] + 2*r*temp$seG[1]*temp$seG[4] + 
                                 2*r*temp$seG[2]*temp$seG[3] + 2*r*temp$seG[2]*temp$seG[4] + 2*r*temp$seG[3]*temp$seG[4])
  } else if(n_exp==3) {
    metadata[metadata$comb_ef == es_markers[i+1],]$meanVarG <- (1/n_exp)^2 * (sum(temp$seG^2, na.rm = TRUE) + 
                                 2*r*temp$seG[1]*temp$seG[2] + 2*r*temp$seG[1]*temp$seG[3]  + 2*r*temp$seG[2]*temp$seG[3])
  } else if (n_exp == 2){
    metadata[metadata$comb_ef == es_markers[i+1],]$meanVarG <- (1/n_exp)^2 * (sum(temp$seG^2, na.rm = TRUE) + 
                                 2*r*temp$seG[1]*temp$seG[2])
  }
}

metadata$G <- metadata$meanG 
metadata$seG <- sqrt(metadata$meanVarG)
metadata <- subset(metadata, select = -c(meanG,meanVarG))

# save data

# write_xlsx(metadata, paste0(directory, "/ES_computed/ES_Computed_Hypothesis.xlsx"))
# write_xlsx(metadata, paste0(directory, "/ES_computed/ES_Computed_Outcome.xlsx"))
write_xlsx(metadata, paste0(directory, "/ES_computed/ES_Computed_nov21.xlsx"))
