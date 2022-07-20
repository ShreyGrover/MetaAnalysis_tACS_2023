##--------------------------------------------------------------------------------------------------------
library("metafor")
library("robumeta")
library("readxl")
library("plotly")

##--------------------------------------------------------------------------------------------------------
directory = "/Users/shreygrover/Desktop/MetaAnalysis_SciTranslMed/Sci_Transl_Med_Submission_2/Analysis/Scripts_Revised/Meta_Datafiles_2021"
mlm.data <- read_excel(paste0(directory, "/Data_for_analysis/Analysis_Outcome_Revised.xlsx"))

# test w/o clinical studies
mlm.data <- mlm.data[mlm.data$population != "Clinical",]

# set to 1 to exclude outliers from the analysis, 0 otherwise
out_exclude = 1

# intensity as a dichotomous variable
mlm.data$intensity_cat <- factor(ifelse(mlm.data$intensity>1, "high", "low"))

# age as a dichotomous variable
mlm.data$age <- factor(ifelse(mlm.data$mean_age>60, "1", "0"))

# creating design variable
mlm.data$design <- factor(ifelse(mlm.data$n == mlm.data$n1, "within", "between"))

# squaring intensity for a post hoc quadratic fit analysis
mlm.data$intensity <- as.numeric(mlm.data$intensity)
mlm.data$intensity_sqr <- mlm.data$intensity^2

# performance and reaction time subsets of data 
mlm.perf <- mlm.data[mlm.data$dv == 'Performance',]
mlm.rt <- mlm.data[mlm.data$dv == 'RT',]

##--------------------------------------------------------------------------------------------------------
# Exclude outliers if needed
if (out_exclude == 1){
  # all
  model_intercept <- robu(formula = y ~ 1, data = mlm.data, studynum = ID, var.eff.size = v, rho = .8, small=FALSE)
  model_rma <- rma.uni(y, v, weights = model_intercept$data.full$r.weights, data=mlm.data)
  rstud <- rstudent(model_rma)
  outliers_all <- model_intercept$data$experiment[abs(rstud$z) >= 1.96]
  mlm.data <- mlm.data[! mlm.data$experiment %in% outliers_all,]
  
  # performance
  model_intercept <- robu(formula = y ~ 1, data = mlm.perf, studynum = ID, var.eff.size = v, rho = .8, small=FALSE)
  model_rma <- rma.uni(y, v, weights = model_intercept$data.full$r.weights, data=mlm.perf)
  rstud <- rstudent(model_rma)
  outliers_perf <- model_intercept$data$experiment[abs(rstud$z) >= 1.96]
  mlm.perf <- mlm.perf[! mlm.perf$experiment %in% outliers_perf,]
  
  # RT
  model_intercept <- robu(formula = y ~ 1, data = mlm.rt, studynum = ID, var.eff.size = v, rho = .8, small=TRUE)
  model_rma <- rma.uni(y, v, weights = model_intercept$data.full$r.weights, data=mlm.rt)
  rstud <- rstudent(model_rma)
  outliers_rt <- model_intercept$data$experiment[abs(rstud$z) >= 1.96]
  mlm.rt <- mlm.rt[! mlm.rt$experiment %in% outliers_rt,]
}

# Run simple meta-regressions for each predictor

# ALL 
## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- robu(formula = y ~ 1 + as.factor(IF), data=mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.numeric(intensity), data= subset(mlm.data, intensity!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + intensity_cat, data= subset(mlm.data, intensity_cat!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(duration), data= mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(online), data= mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(HD), data = mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(current_modeling), data= mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(exploratory), data= mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(design), data= mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(age), data = mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(rs_task), data = mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(neuroguided), data = mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.factor(blinding), data = mlm.data, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_all_output <- rbind(reg_all_output, robu(formula = y ~ 1 + as.numeric(current_density), data = subset(mlm.data, current_density!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)


# PERFORMANCE 

## -------------------------------------------------------------------------------------------------------------------
reg_perf_output <- robu(formula = y ~ 1 + as.factor(IF), data=mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.numeric(intensity), data= subset(mlm.perf, intensity!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + intensity_cat, data= subset(mlm.perf, intensity_cat!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(duration), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(online), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(HD), data = mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(current_modeling), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(exploratory), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(design), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(age), data = mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(rs_task), data = mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(neuroguided), data = mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.factor(blinding), data = mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_perf_output <- rbind(reg_perf_output, robu(formula = y ~ 1 + as.numeric(current_density), data = subset(mlm.perf, current_density!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

# RT

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- robu(formula = y ~ 1 + as.factor(IF), data=mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + intensity, data= subset(mlm.rt, intensity!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## Quadratic fit of intensity to RT effects -----------------------------------------------------------------------

quad <- robu(formula = y ~ 1 + intensity + intensity_sqr, data= subset(mlm.rt, intensity!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + intensity_cat, data= subset(mlm.rt, intensity_cat!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(online), data= mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(HD), data = mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(current_modeling), data= mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(exploratory), data= mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(design), data= mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(rs_task), data = mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(neuroguided), data = mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.factor(blinding), data = mlm.rt, studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------
reg_rt_output <- rbind(reg_rt_output, robu(formula = y ~ 1 + as.numeric(current_density), data = subset(mlm.rt, current_density!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

## ----------------------------------------------------------------------------------------------------------------------------

# save as csv
if (out_exclude==0) {
    write.csv(reg_all_output, paste0(directory, "/Results/WithoutClinical/regression_all.csv"))
    write.csv(reg_perf_output, paste0(directory, "/Results/WithoutClinical/regression_perf.csv"))
    write.csv(reg_rt_output, paste0(directory, "/Results/WithoutClinical/regression_rt.csv"))
}

if (out_exclude == 1){
    write.csv(reg_all_output, paste0(directory, "/Results/WithoutClinical/regression_all_outremoved.csv"))
    write.csv(reg_perf_output, paste0(directory, "/Results/WithoutClinical/regression_perf_outremoved.csv"))
    write.csv(reg_rt_output, paste0(directory, "/Results/WithoutClinical/regression_rt_outremoved.csv"))
}

##--------------------------------------------------------------------------------------------------------

# Interactions 
# During the regression analyses above, we find a significant influence of assessment timing (the "online" variable) on performance effect sizes.
# Accordingly, we now perform subsequent analyses to examine whether the assessment timing variable interacts significantly with any other variable, 
# as long as all interaction combinations of a given pair of variables have at least 10 effect sizes  

# check if enough effects for interaction assessment
table(mlm.perf$online, mlm.perf$IF) # not enough effects for assessing interaction 
table(mlm.perf$online, mlm.perf$duration) 
table(mlm.perf$online, mlm.perf$intensity_cat)
table(mlm.perf$online, mlm.perf$HD)
table(mlm.perf$online, mlm.perf$current_modeling)
table(mlm.perf$online, mlm.perf$age) # not enough effects for assessing interaction 
table(mlm.perf$online, mlm.perf$rs_task) # not enough effects for assessing interaction 
table(mlm.perf$online, mlm.perf$neuroguided) # not enough effects for assessing interaction 

# Run interactions
reg_inter_online <- robu(formula = y ~ 1 + as.factor(online)*intensity, data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table
reg_inter_online <- rbind(reg_inter_online,robu(formula = y ~ 1 + as.factor(online)*as.factor(duration), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)
reg_inter_online <- rbind(reg_inter_online,robu(formula = y ~ 1 + as.factor(online)*as.factor(intensity_cat), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)
reg_inter_online <- rbind(reg_inter_online,robu(formula = y ~ 1 + as.factor(online)*as.factor(HD), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)
reg_inter_online <- rbind(reg_inter_online,robu(formula = y ~ 1 + as.factor(online)*as.factor(current_modeling), data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)
reg_inter_online <- rbind(reg_inter_online,robu(formula = y ~ 1 + as.factor(online)*as.numeric(current_density), data= subset(mlm.perf, current_density!="-"), studynum = ID, v, rho = .8, small = TRUE)$reg_table)

# save as csv
print(reg_inter_online)
if (out_exclude == 0){
  write.csv(reg_inter_online, paste0(directory,"/Results/WithoutClinical/regression_interaction_perf.csv"))
}

if (out_exclude == 1){
  write.csv(reg_inter_online, paste0(directory,"/Results/WithoutClinical/regression_interaction_perf_outremoved.csv"))
}

##--------------------------------------------------------------------------------------------------------

# Probing significance of interactions of categorical variables

# Likelihood Ratio Tests (LRT) (using metafor function rma() that is compatible with R function anova() that performs the LRT)
# to correct for multiplicity we need to divide alpha by 4 (number of interactions we're estimating), thus our alpha is 0.0125

# first fit RVE model to get the weights to use in the subsequent models computed by metafor function rma()
model_intercept <- robu(formula = y ~ 1, data= mlm.perf, studynum = ID, v, rho = .8, small = FALSE)

# fit 2 models: with and without interaction
mod1 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) + as.factor(duration), data = mlm.perf, method="ML")
mod2 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) * as.factor(duration), data = mlm.perf, method="ML")

# compute LRT using anova() function
anova(mod1, mod2)

#
mod1 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) + as.factor(intensity_cat), data = mlm.perf, method="ML")
mod2 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) * as.factor(intensity_cat), data = mlm.perf, method="ML")
anova(mod1, mod2)
#
mod1 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) + as.factor(HD), data = mlm.perf, method="ML")
mod2 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) * as.factor(HD), data = mlm.perf, method="ML")
anova(mod1, mod2)
#
mod1 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) + as.factor(current_modeling), data = mlm.perf, method="ML")
mod2 <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online) * as.factor(current_modeling), data = mlm.perf, method="ML")
anova(mod1, mod2)

##--------------------------------------------------------------------------------------------------------

# Pairwise comparisons

# Current Modeling 
# to correct for multiplicity we need to divide alpha by 2 (number of comparisons we're performing), thus our alpha is 0.025

# this table gives estimates for each effect NOT relative to the intercept, but relative to 0 by using a different syntax (y ~1 + a:b -1) 
reg_inter_pairs <- robu(formula = y ~ 1 + as.factor(online):as.factor(current_modeling)-1, data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table

# use metafor package to perform pairwise comparisons
res <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online):as.factor(current_modeling)-1, data = mlm.perf)

# offline + no modeling vs offline + modeling
anova(res, L=c(1, 0, -1, 0))
# online + no modeling vs online + modeling
anova(res, L=c(0, 1, 0, -1))

# HD (marginally significant LRT)
# to correct for multiplicity we need to divide alpha by 2 (number of comparisons we're performing), thus our alpha is 0.025

reg_inter_pairs <- rbind(reg_inter_pairs, robu(formula = y ~ 1 + as.factor(online):as.factor(HD)-1, data= mlm.perf, studynum = ID, v, rho = .8, small = TRUE)$reg_table)
res <- rma(y, v, weights = model_intercept$data.full$r.weights, mods = ~ as.factor(online):as.factor(HD)-1, data = mlm.perf)

# offline + traditional vs offline + HD
anova(res, L=c(1, 0, -1, 0))
# online + traditional vs online + HD
anova(res, L=c(0, 1, 0, -1))

# save as csv
print(reg_inter_pairs)
if (out_exclude == 0){
  write.csv(reg_inter_pairs, paste0(directory,"/Results/WithoutClinical/regression_interaction_pairwise_perf.csv"))
}

if (out_exclude == 1){
  write.csv(reg_inter_pairs, paste0(directory,"/Results/WithoutClinical/regression_interaction_pairwise_perf_outremoved.csv"))
}

##--------------------------------------------------------------------------------------------------------
