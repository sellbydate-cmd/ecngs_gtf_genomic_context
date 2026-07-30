# Set working directory
setwd("~/R_repo")

install.packages(c("readxl", "ggplot2", "glmmTMB", "emmeans") #installs packages 
                  ,
                 dependencies = TRUE,
                 repos = "https://cloud.r-project.org")
library(readxl)
library(ggplot2)
library(glmmTMB)
library(emmeans) #loads packages

read_basesub <- read_excel("data_for_stats_long.xlsx")
read_regional <- read_excel("regional mutation data.xlsx")
read_muttype <- read_excel("mutation type data.xlsx")
#read data for analyses

read_muttype$mutation_subtype <- relevel(factor(read_muttype$mutation_subtype), ref = "snv") #default reference was complex and that was ultra low frequency, changed to SNV
read_basesub$dose_c <- read_basesub$dose - 1000 #adjusted between -1000-0 to center model effects (-1000 makes dose_c = 0 reflect 1000 mg/kg)

summary(read_basesub$mutation_count)
var(read_basesub$mutation_count)
ggplot(read_basesub, aes(x=mutation_count))+ geom_histogram(color ="white", fill="#4C78A8"+ theme_classic()) 
summary(read_regional$mutation_count)
var(read_regional$mutation_count)
ggplot(read_regional, aes(x=mutation_count))+ geom_histogram(color ="white", fill="#4C78A8"+ theme_classic()) 
summary(read_muttype$mutation_count)
var(read_muttype$mutation_count)
ggplot(read_muttype, aes(x=mutation_count))+ geom_histogram(color ="white", fill="#4C78A8"+ theme_classic()) 
#checking the variance exceeds the mean and using histogram to assess distribution/inform selection of zero-inflated alternatives
#functions could be used to reduce redundancy across this code, not currently implemented


basesub_nb <- glmmTMB(
  mutation_count ~ dose_c * tissue * base_substitution + #*  represents interaction terms + are non-interaction
    offset(log_total_assessed_bases) + 
    (1 | sample), #random effect
  ziformula = ~0, #disables zero inflation, uses NB
  family = nbinom2,
  data = read_basesub,
  control = glmmTMBControl(
    optCtrl = list(iter.max =  10000, eval.max = 10000) #adjusted due to eval limit reached without convergence  
  )
)
basesub_nb$sdr$pdHess #used to check hessian errors

summary(basesub_nb) #gets model outputs
confint(basesub_nb, level = 0.95) #gets confidence intervals
#basesubstitution model

#regional models, non-zero inflated, and zero inflated respectively 
#regional mutation count model
regional_nb <- glmmTMB(
  mutation_count ~ dose + region + tissue + #regional 3 way and 2 way interactions did not work -> NaN (threw non-positive-definite hessian matrix errors), thus simplified 
    offset(log_assessed_bases) + #log offset 
    (1 | sample), #random effect
  ziformula = ~0,
  family = nbinom2,
  data = read_regional
)
regional_nb$sdr$pdHess #used to check hessian errors
summary(regional_nb) #gets model outputs 
confint(regional_nb, level = 0.95) #gets confidence intervals

#regional mutation count model
regional_zinb <- glmmTMB(
  mutation_count ~ dose + region + tissue + #zero inflated model for regional, just did not work (NaN)
    offset(log_assessed_bases) + #log offset 
    (1 | sample), #random effect
  ziformula = ~1, #zero-inflated
  family = nbinom2,
  data = read_regional #specifies dataset for model
)

summary(regional_zinb) #gets model output
regional_zinb$sdr$pdHess #used to check hess errors
#NOTE. Histograms indicated abundance of 0s, thus the use of zero inflated were attempted and AICs were compared
AIC(regional_nb, regional_zinb) #gets AICs to compare models



#code for mutation subtype comparisons (zero inflated, non-zero inflated neg binomial, respectively)
#NOTE. Histograms indicated abundance of 0s, thus the use of zero inflated were attempted and AICs were compared
muttype_nb <- glmmTMB(
  mutation_count ~ dose + tissue * mutation_subtype + #full interaction did not work properly, had to simplify (Nb model performed best, based of AIC)
    offset(log_total_assessed_bases) + #log offset 
    (1 | sample), #random effect
    ziformula = ~0, #disables zero-inflation
  family = nbinom2,
  data = read_muttype
)

summary(muttype_nb) #model output
confint(muttype_nb, level = 0.95) #for CIs 
muttype_nb$sdr$pdHess #hess errors check

muttype_zinb <- glmmTMB(
  mutation_count ~ dose + tissue * mutation_subtype +
    offset(log_total_assessed_bases) + #log offset 
    (1 | sample), #random effect
  ziformula = ~1, #zero inflated
  family = nbinom2,
  data = read_muttype
)
summary(muttype_zinb)
confint(muttype_zinb)
muttype_zinb$sdr$pdHess #hess error check
muttype_zinb$fit$convergence == 0 #checks convergence

#summary
AIC(muttype_nb, muttype_zinb) #used to compare AICs between models

#post-hoc testing

emm_region <- emmeans(regional_nb, ~ region) #could not check between tissues since interaction was not included
pairs(emm_region, adjust = "tukey", type = "response", infer = c(TRUE, TRUE)) #response (transforms from log scale by taking exponentiated log estimate)
#also infer = 95% CI by default

#remove dose_c centering higher in the document to ensure this reflects baseline response
emm_base_sub <- emmeans(basesub_nb, ~ base_substitution | tissue, at = list(dose_c = 0)) #note, check assigned dose_c value since this uses dose_c and it is set higher in document
pairs(emm_base_sub, adjust = "tukey", type = "response", infer = c(TRUE, TRUE)) 

#remove dose_c centering higher in the document to ensure this reflects maximum response
emm_base_sub_1000 <- emmeans(basesub_nb, ~ base_substitution | tissue, at = list(dose_c = 1000))
pairs(emm_base_sub_1000, adjust = "tukey", type = "response", infer = c(TRUE, TRUE))

trends <- emtrends(basesub_nb, ~ tissue * base_substitution, var = "dose_c")
summary(pairs(trends, by = "base_substitution", reverse = TRUE, infer = c(TRUE, TRUE)), by = NULL, adjust= "bonferroni")
#reverse = true (makes liver the reference), according to emmeans documentation no adjustment is applied if comparing two means

emm_muttype <- emmeans(muttype_nb, ~ mutation_subtype | tissue)
pairs(emm_muttype, adjust = "tukey", type = "response", infer = c(TRUE, TRUE))

citation("glmmTMB") 
citation("readxl")
citation("ggplot2")
citation("emmeans")
#gets citations for packages
