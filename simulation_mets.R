# (first run script "parameterization")
# load libraries ----------------------------------------------------------
library(MicSim)
library(tidyverse)
library(writexl)

# setting path to data files
file_path <- "./data_files/" 

# defining initial pop MetS -----------------------------------------------
# setting sample size to 500K
N <- 10^5 * 5

# creating birth cohort 2000 (they enter simulation in 2018 at age 18)
set.seed(20210607)
initBirthDatesRange <- chron(dates=c("1/1/2000","31/12/2000"), format=c(dates="d/m/Y"),out.format=c(dates="d/m/year"))
birthDates <- dates(initBirthDatesRange[1] + runif(N, min=0, max=diff(initBirthDatesRange)))

# setting minimum and maximum age 
minage <- 18
maxage <- 65

# setting simulation horizon
simHorizon <- setSimHorizon(startDate="01/01/2018", endDate="31/12/2065")

# initial states: intersection between gender, education level and health status
states_mets <- c("noms.lo.m","noms.lo.f","noms.hi.m","noms.hi.f",
                 "ms.lo.m","ms.lo.f","ms.hi.m","ms.hi.f")

# absorbing states: dead (must be included even though mortality rate set to zero)
absStates <- "dead"

# all states 
stateSpace_mets <- expand.grid(states_mets)

# prevalence of MetS in starting pop depends on a stochastic process driven by 
# percentages provided by regressions. assume equal distribution among men and 
# women and high and low education
m.lo.mets <- sample(x=c("ms.lo.m","noms.lo.m"),
               prob=c(.0238,1-.0238),size=N/4,replace=TRUE) 

f.lo.mets <- sample(x=c("ms.lo.f","noms.lo.f"),
               prob=c(.0200,1-.0200),size=N/4,replace=TRUE)

m.hi.mets <- sample(x=c("ms.hi.m","noms.hi.m"),
               prob=c(.0068,1-.0068),size=N/4,replace=TRUE)

f.hi.mets <- sample(x=c("ms.hi.f","noms.hi.f"),
               prob=c(.0057,1-.0057),size=N/4,replace=TRUE)

initStates_mets <- c(m.lo.mets, f.lo.mets, m.hi.mets, f.hi.mets)

initPop_mets <- data.frame(ID=1:N, birthDate=birthDates,initState=initStates_mets)


# absorbing state for all models. this is needed even though we don't account for mortality
absTransitions <- c("dead","mortRates")


# creating transition matrices --------------------------------------------
# MetS observed data: transition pattern and assignment of functions specifying transition rates
healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                               "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                               "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                               "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                             c("noms_ms.lo.m","noms_ms.lo.f",
                               "noms_ms.hi.m","noms_ms.hi.f",
                               "ms_noms.lo.m","ms_noms.lo.f",
                               "ms_noms.hi.m","ms_noms.hi.f"))
allTransitions_mets <- healthTrMatrix_mets
transitionMatrix_mets <- buildTransitionMatrix(allTransitions=allTransitions_mets,
                                               stateSpace=stateSpace_mets, 
                                               absTransitions=absTransitions)

# MetS counterfactual data smk: transition pattern and assignment of functions specifying transition rates
cf_smk_healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                                      "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                                      "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                                      "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                                    c("cf_smk_noms_ms.lo.m","cf_smk_noms_ms.lo.f",
                                      "cf_smk_noms_ms.hi.m","cf_smk_noms_ms.hi.f",
                                      "cf_smk_ms_noms.lo.m","cf_smk_ms_noms.lo.f",
                                      "cf_smk_ms_noms.hi.m","cf_smk_ms_noms.hi.f"))
cf_smk_allTransitions_mets <- cf_smk_healthTrMatrix_mets
cf_smk_transitionMatrix_mets <- buildTransitionMatrix(allTransitions=cf_smk_allTransitions_mets,
                                                      stateSpace=stateSpace_mets,
                                                      absTransitions=absTransitions)

# MetS counterfactual data alc: transition pattern and assignment of functions specifying transition rates
cf_alc_healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                                      "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                                      "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                                      "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                                    c("cf_alc_noms_ms.lo.m","cf_alc_noms_ms.lo.f",
                                      "cf_alc_noms_ms.hi.m","cf_alc_noms_ms.hi.f",
                                      "cf_alc_ms_noms.lo.m","cf_alc_ms_noms.lo.f",
                                      "cf_alc_ms_noms.hi.m","cf_alc_ms_noms.hi.f"))
cf_alc_allTransitions_mets <- cf_alc_healthTrMatrix_mets
cf_alc_transitionMatrix_mets <- buildTransitionMatrix(allTransitions=cf_alc_allTransitions_mets,
                                                      stateSpace=stateSpace_mets,
                                                      absTransitions=absTransitions)

# MetS counterfactual data diet: transition pattern and assignment of functions specifying transition rates
cf_diet_healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                                       "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                                       "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                                       "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                                     c("cf_diet_noms_ms.lo.m","cf_diet_noms_ms.lo.f",
                                       "cf_diet_noms_ms.hi.m","cf_diet_noms_ms.hi.f",
                                       "cf_diet_ms_noms.lo.m","cf_diet_ms_noms.lo.f",
                                       "cf_diet_ms_noms.hi.m","cf_diet_ms_noms.hi.f"))
cf_diet_allTransitions_mets <- cf_diet_healthTrMatrix_mets
cf_diet_transitionMatrix_mets <- buildTransitionMatrix(allTransitions=cf_diet_allTransitions_mets,
                                                       stateSpace=stateSpace_mets,
                                                       absTransitions=absTransitions)

# MetS counterfactual data hl: transition pattern and assignment of functions specifying transition rates
cf_hl_healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                                     "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                                     "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                                     "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                                   c("cf_hl_noms_ms.lo.m","cf_hl_noms_ms.lo.f",
                                     "cf_hl_noms_ms.hi.m","cf_hl_noms_ms.hi.f",
                                     "cf_hl_ms_noms.lo.m","cf_hl_ms_noms.lo.f",
                                     "cf_hl_ms_noms.hi.m","cf_hl_ms_noms.hi.f"))
cf_hl_allTransitions_mets <- cf_hl_healthTrMatrix_mets
cf_hl_transitionMatrix_mets <- buildTransitionMatrix(allTransitions=cf_hl_allTransitions_mets,
                                                     stateSpace=stateSpace_mets,
                                                     absTransitions=absTransitions)

# MetS counterfactual data joint: transition pattern and assignment of functions specifying transition rates
cf_joint_healthTrMatrix_mets <- cbind(c("noms.lo.m->ms.lo.m","noms.lo.f->ms.lo.f",
                                        "noms.hi.m->ms.hi.m","noms.hi.f->ms.hi.f",
                                        "ms.lo.m->noms.lo.m","ms.lo.f->noms.lo.f",
                                        "ms.hi.m->noms.hi.m","ms.hi.f->noms.hi.f"),
                                      c("cf_joint_noms_ms.lo.m","cf_joint_noms_ms.lo.f",
                                        "cf_joint_noms_ms.hi.m","cf_joint_noms_ms.hi.f",
                                        "cf_joint_ms_noms.lo.m","cf_joint_ms_noms.lo.f",
                                        "cf_joint_ms_noms.hi.m","cf_joint_ms_noms.hi.f"))
cf_joint_allTransitions_mets <- cf_joint_healthTrMatrix_mets
cf_joint_transitionMatrix_mets <- buildTransitionMatrix(allTransitions=cf_joint_allTransitions_mets,
                                                        stateSpace=stateSpace_mets,
                                                        absTransitions=absTransitions)


# running simulations -----------------------------------------------------
# supplying the micsim function with the required input created in the code above
simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=transitionMatrix_mets, 
                              absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                              cores = 80, seeds = 20210607) 

cf_smk_simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=cf_smk_transitionMatrix_mets, 
                                     absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                                     cores = 80, seeds = 20210607) 

cf_alc_simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=cf_alc_transitionMatrix_mets, 
                                     absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                                     cores = 80, seeds = 20210607) 

cf_diet_simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=cf_diet_transitionMatrix_mets, 
                                      absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                                      cores = 80, seeds = 20210607) 

cf_hl_simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=cf_hl_transitionMatrix_mets, 
                                    absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                                    cores = 80, seeds = 20210607) 

cf_joint_simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=cf_joint_transitionMatrix_mets, 
                                       absStates=absStates, maxAge=maxage, simHorizon=simHorizon, 
                                       cores = 80, seeds = 20210607)

# converting output to spell format ---------------------------------------
# preparing datasets: prefix "cf_XXX_" gives which counterfactual this is

# saving/loading output of simulation, so it doesn't have to be rerun
save(simpop_mets, file = paste0(file_path, "simpop_mets.rda"))
save(cf_smk_simpop_mets, file = paste0(file_path, "cf_smk_simpop_mets.rda"))
save(cf_alc_simpop_mets, file = paste0(file_path, "cf_alc_simpop_mets.rda"))
save(cf_diet_simpop_mets, file = paste0(file_path, "cf_diet_simpop_mets.rda"))
save(cf_hl_simpop_mets, file = paste0(file_path, "cf_hl_simpop_mets.rda"))
save(cf_joint_simpop_mets, file = paste0(file_path, "cf_joint_simpop_mets.rda"))

load(paste0(file_path, "simpop_mets.rda"))
load(paste0(file_path, "cf_smk_simpop_mets.rda"))
load(paste0(file_path, "cf_alc_simpop_mets.rda"))
load(paste0(file_path, "cf_diet_simpop_mets.rda"))
load(paste0(file_path, "cf_hl_simpop_mets.rda"))
load(paste0(file_path, "cf_joint_simpop_mets.rda"))

# first duplicating data
msm_mets <- simpop_mets %>% arrange(ID, transitionTime)
cf_smk_msm_mets <- cf_smk_simpop_mets %>% arrange(ID, transitionTime)
cf_alc_msm_mets <- cf_alc_simpop_mets %>% arrange(ID, transitionTime)
cf_diet_msm_mets <- cf_diet_simpop_mets %>% arrange(ID, transitionTime)
cf_hl_msm_mets <- cf_hl_simpop_mets %>% arrange(ID, transitionTime)
cf_joint_msm_mets <- cf_joint_simpop_mets %>% arrange(ID, transitionTime)

# creating index to identify each transition
msm_mets <- msm_mets %>% group_by(ID) %>% mutate(index=1:n()) %>% ungroup()
cf_smk_msm_mets <- cf_smk_msm_mets %>% group_by(ID) %>% mutate(index=1:n())  %>% ungroup()
cf_alc_msm_mets <- cf_alc_msm_mets %>% group_by(ID) %>% mutate(index=1:n())  %>% ungroup()
cf_diet_msm_mets <- cf_diet_msm_mets %>% group_by(ID) %>% mutate(index=1:n())  %>% ungroup()
cf_hl_msm_mets <- cf_hl_msm_mets %>% group_by(ID) %>% mutate(index=1:n())  %>% ungroup()
cf_joint_msm_mets <- cf_joint_msm_mets %>% group_by(ID) %>% mutate(index=1:n())  %>% ungroup()

# only selecting first row per individual and adding labels for 
# gender, education and health status at start of model
obs1_mets <- msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

cf_smk_obs1_mets <- cf_smk_msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

cf_alc_obs1_mets <- cf_alc_msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

cf_diet_obs1_mets <- cf_diet_msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

cf_hl_obs1_mets <- cf_hl_msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

cf_joint_obs1_mets <- cf_joint_msm_mets %>%
  filter(index == 1) %>%
  mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
         education = if_else(str_detect(initState, ".lo."), "Low", "High"),
         health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
         agestart = minage)

# recoding health status 
msm_mets$health <-  if_else(str_detect(msm_mets$To, "noms."), "NoMS", "MS")
cf_smk_msm_mets$health <-  if_else(str_detect(cf_smk_msm_mets$To, "noms."), "NoMS", "MS")
cf_alc_msm_mets$health <-  if_else(str_detect(cf_alc_msm_mets$To, "noms."), "NoMS", "MS")
cf_diet_msm_mets$health <-  if_else(str_detect(cf_diet_msm_mets$To, "noms."), "NoMS", "MS")
cf_hl_msm_mets$health <-  if_else(str_detect(cf_hl_msm_mets$To, "noms."), "NoMS", "MS")
cf_joint_msm_mets$health <-  if_else(str_detect(cf_joint_msm_mets$To, "noms."), "NoMS", "MS")

# rearrange data from msm format to spell format (one line per state)
# mets observed
obs1_mets <- obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

gendereduc_mets <- select(obs1_mets,ID,gender,education)

obs1_mets <- obs1_mets %>% 
  select(ID,agestart,agestop,health)

otherobs_mets <- list()
for (i in 1:max(msm_mets$index)){
  a <- select(filter(msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

otherobs_mets <- bind_rows(otherobs_mets)
otherobs_mets$agestop <- ifelse(is.na(otherobs_mets$agestop),maxage,otherobs_mets$agestop)

spell_mets <- bind_rows(obs1_mets,otherobs_mets)
spell_mets <- spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# mets counterfactual smk
cf_smk_obs1_mets <- cf_smk_obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

cf_smk_gendereduc_mets <- select(cf_smk_obs1_mets,ID,gender,education)

cf_smk_obs1_mets <- select(cf_smk_obs1_mets,ID,agestart,agestop,health)

cf_smk_otherobs_mets <- list()
for (i in 1:max(cf_smk_msm_mets$index)){
  a <- select(filter(cf_smk_msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(cf_smk_msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  cf_smk_otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

cf_smk_otherobs_mets <- bind_rows(cf_smk_otherobs_mets)
cf_smk_otherobs_mets$agestop <- ifelse(is.na(cf_smk_otherobs_mets$agestop),maxage,cf_smk_otherobs_mets$agestop)

cf_smk_spell_mets <- bind_rows(cf_smk_obs1_mets,cf_smk_otherobs_mets)
cf_smk_spell_mets <- cf_smk_spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# mets counterfactual alc
cf_alc_obs1_mets <- cf_alc_obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

cf_alc_gendereduc_mets <- select(cf_alc_obs1_mets,ID,gender,education)

cf_alc_obs1_mets <- select(cf_alc_obs1_mets,ID,agestart,agestop,health)

cf_alc_otherobs_mets <- list()
for (i in 1:max(cf_alc_msm_mets$index)){
  a <- select(filter(cf_alc_msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(cf_alc_msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  cf_alc_otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

cf_alc_otherobs_mets <- bind_rows(cf_alc_otherobs_mets)
cf_alc_otherobs_mets$agestop <- ifelse(is.na(cf_alc_otherobs_mets$agestop),maxage,cf_alc_otherobs_mets$agestop)

cf_alc_spell_mets <- bind_rows(cf_alc_obs1_mets,cf_alc_otherobs_mets)
cf_alc_spell_mets <- cf_alc_spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# mets counterfactual diet
cf_diet_obs1_mets <- cf_diet_obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

cf_diet_gendereduc_mets <- select(cf_diet_obs1_mets,ID,gender,education)

cf_diet_obs1_mets <- select(cf_diet_obs1_mets,ID,agestart,agestop,health)

cf_diet_otherobs_mets <- list()
for (i in 1:max(cf_diet_msm_mets$index)){
  a <- select(filter(cf_diet_msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(cf_diet_msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  cf_diet_otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

cf_diet_otherobs_mets <- bind_rows(cf_diet_otherobs_mets)
cf_diet_otherobs_mets$agestop <- ifelse(is.na(cf_diet_otherobs_mets$agestop),maxage,cf_diet_otherobs_mets$agestop)

cf_diet_spell_mets <- bind_rows(cf_diet_obs1_mets,cf_diet_otherobs_mets)
cf_diet_spell_mets <- cf_diet_spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# mets counterfactual hl
cf_hl_obs1_mets <- cf_hl_obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

cf_hl_gendereduc_mets <- select(cf_hl_obs1_mets,ID,gender,education)

cf_hl_obs1_mets <- select(cf_hl_obs1_mets,ID,agestart,agestop,health)

cf_hl_otherobs_mets <- list()
for (i in 1:max(cf_hl_msm_mets$index)){
  a <- select(filter(cf_hl_msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(cf_hl_msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  cf_hl_otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

cf_hl_otherobs_mets <- bind_rows(cf_hl_otherobs_mets)
cf_hl_otherobs_mets$agestop <- ifelse(is.na(cf_hl_otherobs_mets$agestop),maxage,cf_hl_otherobs_mets$agestop)

cf_hl_spell_mets <- bind_rows(cf_hl_obs1_mets,cf_hl_otherobs_mets)
cf_hl_spell_mets <- cf_hl_spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# mets counterfactual joint
cf_joint_obs1_mets <- cf_joint_obs1_mets %>%
  mutate(agestop = if_else(is.na(transitionAge), maxage, transitionAge))

cf_joint_gendereduc_mets <- select(cf_joint_obs1_mets,ID,gender,education)

cf_joint_obs1_mets <- select(cf_joint_obs1_mets,ID,agestart,agestop,health)

cf_joint_otherobs_mets <- list()
for (i in 1:max(cf_joint_msm_mets$index)){
  a <- select(filter(cf_joint_msm_mets,index==i),ID,transitionAge,health)
  a <- rename(a,agestart=transitionAge)
  b <- select(filter(cf_joint_msm_mets,index==i+1),ID,transitionAge)
  b <- rename(b,agestop=transitionAge)
  
  cf_joint_otherobs_mets[[i]] <- left_join(a,b,by="ID")
}

cf_joint_otherobs_mets <- bind_rows(cf_joint_otherobs_mets)
cf_joint_otherobs_mets$agestop <- ifelse(is.na(cf_joint_otherobs_mets$agestop),maxage,cf_joint_otherobs_mets$agestop)

cf_joint_spell_mets <- bind_rows(cf_joint_obs1_mets,cf_joint_otherobs_mets)
cf_joint_spell_mets <- cf_joint_spell_mets %>%
  left_join(gendereduc_mets, by = "ID") %>%
  arrange(ID, agestart) %>%
  drop_na(.) %>%
  mutate(length = agestop - agestart)

# assign index for each state per individual
spell_mets <- spell_mets %>% group_by(ID) %>% mutate(index=1:n())
cf_smk_spell_mets <- cf_smk_spell_mets %>% group_by(ID) %>% mutate(index=1:n()) 
cf_alc_spell_mets <- cf_alc_spell_mets %>% group_by(ID) %>% mutate(index=1:n())
cf_diet_spell_mets <- cf_diet_spell_mets %>% group_by(ID) %>% mutate(index=1:n())
cf_hl_spell_mets <- cf_hl_spell_mets %>% group_by(ID) %>% mutate(index=1:n())
cf_joint_spell_mets <- cf_joint_spell_mets %>% group_by(ID) %>% mutate(index=1:n())

# merging results from the simulations and creating a variable (sim) to 
# differentiate the simulations. table order variable is meant to help with
# formatting of output when summarizing results later on
mydata_mets <- bind_rows(mutate(spell_mets, sim = "observed"), 
                         mutate(cf_smk_spell_mets, sim = "counterfactual smoking"),
                         mutate(cf_alc_spell_mets, sim = "counterfactual alcohol"),
                         mutate(cf_diet_spell_mets, sim = "counterfactual diet"),
                         mutate(cf_hl_spell_mets, sim = "counterfactual health literacy"),
                         mutate(cf_joint_spell_mets, sim = "counterfactual joint")) %>%
  mutate(sim = factor(sim, 
                levels = c("observed", "counterfactual smoking",
                           "counterfactual alcohol", "counterfactual diet",
                           "counterfactual health literacy", "counterfactual joint")),
         table_order = case_when(sim == "observed" & education == "High" ~ 0,
                           sim == "observed" & education == "Low"  ~ 1,
                           sim == "counterfactual smoking"  & education == "Low"  ~ 2,
                           sim == "counterfactual alcohol"  & education == "Low"  ~ 3,
                           sim == "counterfactual diet"  & education == "Low"  ~ 4,
                           sim == "counterfactual health literacy"  & education == "Low"  ~ 5,
                           sim == "counterfactual joint"  & education == "Low"  ~ 6))

# creating dataset for sensitivty analyses. Must stay in a state for at the least
# 6 months for it to be valid transition. 
mydata_mets_sensitivity <- mydata_mets %>%
  group_by(sim, ID) %>%
  mutate(
    health = case_when(
      #accounting for three consecutive transitions < 6 months
      health == "NoMS" & index == 1 & length < 0.5 & 
        lead(length) < 0.5 & lead(length, n = 2) < 0.5 ~ "MS",
      health == "MS" & index == 1 & length < 0.5 & 
        lead(length) < 0.5 & lead(length, n = 2) < 0.5 ~ "NoMS",
      index == 2 & length < 0.5 & lag(length) < 0.5 &
        lead(length) < 0.5 ~ health,
      length < 0.5 & lag(length) < 0.5 & lead(length) < 0.5 ~ health,
      health == "NoMS" & length < 0.5 & lag(length) < 0.5 &
        lag(length, n = 2) < 0.5 ~ "MS",
      health == "MS" & length < 0.5 & lag(length) < 0.5 &
        lag(length, n = 2) < 0.5 ~ "NoMS",
      #accounting for two consecutive transitions < 6 months
      index == 1 & length < 0.5 & lead(length) < 0.5 ~ health,
      health == "NoMS" & index == 2 & length < 0.5 & 
        lag(length) < 0.5 ~ "MS",
      health == "MS" & index == 2 & length < 0.5 & 
        lag(length) < 0.5 ~ "NoMS",
      # line below also accounts for middle transition when there are 3 consecutive transitions
      length < 0.5 & lag(length) < 0.5 ~ health, 
      health == "MS" & length < 0.5 & 
        lead(length) < 0.5 ~ "NoMS", 
      health == "NoMS" & length < 0.5 & 
        lead(length) < 0.5 ~ "MS", 
      health == "MS"   & length < 0.5 ~ "NoMS",
      health == "NoMS" & length < 0.5 ~ "MS",
      T ~ health),
    unique_transition = case_when(
      health != lead(health) & !(is.na(lead(health))) &
        health != lag(health)  & !(is.na(lag(health))) ~ 1,
      health == lag(health)  & is.na(lead(health))     ~ 1,
      health != lag(health)  & is.na(lead(health))     ~ 1,
      health != lead(health) & !(is.na(lead(health)))  ~ 1,
      n() == 1                                         ~ 1,
      T ~ 0)) %>%
  add_tally(unique_transition, name = "n_transitions") %>%
  ungroup()

# exporting data so I can reload it later on
save(mydata_mets, file = paste0(file_path, "mydata_mets.rda"))
save(mydata_mets_sensitivity, file = paste0(file_path, "mydata_mets_sensitivity.rda"))

load(paste0(file_path, "mydata_mets.rda"))
load(paste0(file_path, "mydata_mets_sensitivity.rda"))
mydata_mets <- mydata_mets_fixed
# analyze results ---------------------------------------------------------
# prevalence --------------------------------------------------------------
# life course prevalence aka proportion who ever had MetS filtering
# by index <= 2 to avoid duplicate individuals in primary analysis. 
# filtering slightly differently in sensitivity analysis because it 
# is possible for the first instance of MetS to after the second index
lc_prev_mets <- mydata_mets %>%
  filter(index <= 2, health == "MS", !is.na(table_order)) %>% 
  group_by(sim, education, gender) %>%
  summarise(proportion_MetS = n()/(N*.25) * 100, .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (n%)", values_from = proportion_MetS) %>%
  arrange(sim); lc_prev_mets 

lc_prev_mets_sens <- mydata_mets_sensitivity %>%
  filter(health == "MS", !is.na(table_order)) %>% 
  distinct(ID, sim, .keep_all = T) %>%
  group_by(sim, education, gender) %>%
  summarise(proportion_MetS = n()/(N*.25) * 100, .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (n%)", values_from = proportion_MetS) %>%
  arrange(sim); lc_prev_mets_sens

# reformatting data frames to use in bar plots later on
lc_prev_mets_diff <- lc_prev_mets %>%
  mutate(across(where(is.numeric), ~ .x - first(.))) %>%
  filter(education == "Low") %>% 
  rename(Female = `Female (n%)`, Male = `Male (n%)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "lc_prev"); lc_prev_mets_diff 

lc_prev_mets_sens_diff <- lc_prev_mets_sens %>%
  mutate(across(where(is.numeric), ~ .x - first(.))) %>%
  filter(education == "Low") %>%
  rename(Female = `Female (n%)`, Male = `Male (n%)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "lc_prev"); lc_prev_mets_sens_diff 

lc_prev_mets_both <- lc_prev_mets_diff %>%
  mutate(analysis = "Primary") %>%
  bind_rows(mutate(lc_prev_mets_sens_diff, analysis = "Sensitivity")); lc_prev_mets_both

rm(lc_prev_mets_diff, lc_prev_mets_sens_diff)


# plotting the prevalence of MetS at each age by sex and education
# using observed data in primary analysis
tmp_data <- mydata_mets %>%
  filter(sim == "observed")

prev_mets <- list()
for (i in minage:maxage){
  prev_mets[[i]] <- tmp_data %>%
    filter((health == "MS" & agestart <= i & agestop >= i)) %>%
    group_by(sim, gender, education) %>%
    summarise(prev = n()/(N*.25) * 100, .groups = "keep") %>%
    mutate(age = i) %>%
    ungroup()
}

# merging prevalence estimates into one data frame
prev_mets <- bind_rows(prev_mets) 

# plotting the age-specific prevalence of MetS stratified
# by education and sex
ggplot(prev_mets) +
  geom_line(aes(age, prev, color=gender, 
                linetype = education), size = 1)+
  labs(y = "Prevalence of MetS %",
       x = "Age (years)",
       colour = "Sex",
       linetype = "Education") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.box.margin = margin(0,0,0,0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot",
        plot.margin = margin(0,0,0,0))
ggsave("./plots/age_specific_prev_mets.pdf", plot = last_plot(),
       width = 6, height = 4, units = "in")

rm(tmp_data)

# age of onset ------------------------------------------------------------
# mean age at onset if no MS age 18 filtering by index == 2 to 
# avoid duplicate individuals and exclude those who had MetS when they
# entered the cohort. filtering slightly differently in sensitivity analysis 
# first creating initstate variable to identify health status at age 18
# then selecting individuals who did not have MetS at age 18 (when entering cohort),
# but had MetS at some point. Then using distinct to only keep the first occurrence
# of MetS.
mean_onset_mets <- mydata_mets %>%
  filter(index == 2, health == "MS", !is.na(table_order)) %>%
  group_by(sim, education, gender) %>%
  summarise(age_onset = mean(agestart), .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (years)", values_from = age_onset) %>%
  arrange(sim); mean_onset_mets

mean_onset_mets_sens <- mydata_mets_sensitivity %>%
  filter(!is.na(table_order)) %>%
  group_by(sim, ID) %>%
  mutate(initstate = first(health)) %>%
  filter(initstate == "NoMS" & health == "MS") %>%
  distinct(ID, sim, .keep_all = T) %>%
  group_by(sim, education, gender) %>%
  summarise(age_onset = mean(agestart), .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (years)", values_from = age_onset) %>%
  arrange(sim); mean_onset_mets_sens

# reformatting data frames to use in bar plots later on
mean_onset_mets_diff <- mean_onset_mets %>%
  mutate(across(where(is.numeric), ~ (.x - first(.)) * -1)) %>%
  filter(education == "Low") %>% 
  rename(Female = `Female (years)`, Male = `Male (years)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "mean_onset"); mean_onset_mets_diff 

mean_onset_mets_sens_diff <- mean_onset_mets_sens %>%
  mutate(across(where(is.numeric), ~ (.x - first(.)) * -1)) %>%
  filter(education == "Low") %>%
  rename(Female = `Female (years)`, Male = `Male (years)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "mean_onset"); mean_onset_mets_sens_diff 

mean_onset_mets_both <- mean_onset_mets_diff %>%
  mutate(analysis = "Primary") %>%
  bind_rows(mutate(mean_onset_mets_sens_diff, analysis = "Sensitivity"))

rm(mean_onset_mets_diff, mean_onset_mets_sens_diff)


# mean duration of MetS ---------------------------------------------------
# mean duration of MetS between ages 18 and 65
# simply summing up total time spent with MetS and dividing by number of
# unique individuals who ever had MetS in each strata.
duration_of_mets <- mydata_mets %>%
  filter(health == "MS", !is.na(table_order)) %>%
  group_by(sim, education, gender) %>%
  summarise(duration_MetS = sum(length)/length(unique(ID)), .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (years)", values_from = duration_MetS) %>%
  arrange(sim); duration_of_mets

duration_of_mets_sens <- mydata_mets_sensitivity %>%
  filter(health == "MS", !is.na(table_order)) %>%
  group_by(sim, education, gender) %>%
  summarise(duration_MetS = sum(length)/length(unique(ID)), .groups = "drop") %>%
  pivot_wider(names_from = gender, names_glue = "{gender} (years)", values_from = duration_MetS) %>%
  arrange(sim); duration_of_mets_sens

# reformatting data frames to use in bar plots later on
duration_of_mets_diff <- duration_of_mets %>%
  mutate(across(where(is.numeric), ~ .x - first(.))) %>%
  filter(education == "Low") %>% 
  rename(Female = `Female (years)`, Male = `Male (years)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "duration"); duration_of_mets_diff 

duration_of_mets_sens_diff <- duration_of_mets_sens %>%
  mutate(across(where(is.numeric), ~ .x - first(.))) %>%
  filter(education == "Low") %>%
  rename(Female = `Female (years)`, Male = `Male (years)`) %>%
  pivot_longer(cols = c(Female, Male), names_to = "sex", values_to = "duration"); duration_of_mets_sens_diff 

duration_of_mets_both <- duration_of_mets_diff %>%
  mutate(analysis = "Primary") %>%
  bind_rows(mutate(duration_of_mets_sens_diff, analysis = "Sensitivity")) 

rm(duration_of_mets_diff, duration_of_mets_sens_diff) 


# barplots for differences between groups ---------------------------------
# primary analysis
lc_prev_dif_plot <- lc_prev_mets_both %>%
  filter(analysis == "Primary") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                         sim == "counterfactual diet" ~ "Counterfactual diet quality",
                         sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                         sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                         sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                         sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(lc_prev_percent = if_else(sim == "Observed data", "Reference",
                                   sprintf("-%1.1f%%", lc_prev_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = lc_prev, y =  reorder(sim, -lc_prev)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = lc_prev, y =  reorder(sim, -lc_prev),
                label = lc_prev_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in life course prevalence of MetS (% points)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.75)); lc_prev_dif_plot
ggsave("./plots/lc_prev_dif_plot.pdf", plot = lc_prev_dif_plot)

mean_onset_dif_plot <- mean_onset_mets_both %>%
  filter(analysis == "Primary") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                                   sim == "counterfactual diet" ~ "Counterfactual diet quality",
                                   sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                                   sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                                   sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                                   sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(mean_onset_percent = if_else(sim == "Observed data", "Reference",
                                      sprintf("-%1.1f%%", mean_onset_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = mean_onset, y =  reorder(sim, -mean_onset)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = mean_onset, y =  reorder(sim, -mean_onset),
                label = mean_onset_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in mean age of onset of MetS (years)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank()); mean_onset_dif_plot
ggsave("./plots/mean_onset_dif_plot.pdf", plot = mean_onset_dif_plot)

mean_duration_dif_plot <- duration_of_mets_both %>%
  filter(analysis == "Primary") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                                   sim == "counterfactual diet" ~ "Counterfactual diet quality",
                                   sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                                   sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                                   sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                                   sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(duration_percent = if_else(sim == "Observed data", "Reference",
                                    sprintf("-%1.1f%%", duration_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = duration, y =  reorder(sim, -duration)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = duration, y =  reorder(sim, -duration),
                label = duration_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in mean duration of MetS (years)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank()); mean_duration_dif_plot
ggsave("./plots/mean_duration_dif_plot.pdf", plot = mean_duration_dif_plot)


#sensitivity analysis
lc_prev_dif_plot_sens <- lc_prev_mets_both %>%
  filter(analysis == "Sensitivity") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                                   sim == "counterfactual diet" ~ "Counterfactual diet quality",
                                   sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                                   sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                                   sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                                   sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(lc_prev_percent = if_else(sim == "Observed data", "Reference",
                                   sprintf("-%1.1f%%", lc_prev_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = lc_prev, y =  reorder(sim, -lc_prev)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = lc_prev, y =  reorder(sim, -lc_prev),
                label = lc_prev_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in life course prevalence of MetS (% points)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(hjust = 0.75)); lc_prev_dif_plot_sens
ggsave("./plots/lc_prev_dif_plot_sens.pdf", plot = lc_prev_dif_plot_sens)

mean_onset_dif_plot_sens <- mean_onset_mets_both %>%
  filter(analysis == "Sensitivity") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                                   sim == "counterfactual diet" ~ "Counterfactual diet quality",
                                   sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                                   sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                                   sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                                   sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(mean_onset_percent = if_else(sim == "Observed data", "Reference",
                                      sprintf("-%1.1f%%", mean_onset_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = mean_onset, y =  reorder(sim, -mean_onset)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = mean_onset, y =  reorder(sim, -mean_onset),
                label = mean_onset_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in mean age of onset of MetS (years)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank()); mean_onset_dif_plot_sens
ggsave("./plots/mean_onset_dif_plot_sens.pdf", plot = mean_onset_dif_plot_sens)

mean_duration_of_dif_plot_sens <- duration_of_mets_both %>%
  filter(analysis == "Sensitivity") %>%
  mutate(sim = as.factor(case_when(sim == "observed" ~ "Observed data",
                                   sim == "counterfactual diet" ~ "Counterfactual diet quality",
                                   sim == "counterfactual health literacy" ~ "Counterfactual health literacy",
                                   sim == "counterfactual alcohol" ~ "Counterfactual alcohol intake",
                                   sim == "counterfactual smoking" ~ "Counterfactual smoking behaviour",
                                   sim == "counterfactual joint" ~ "Counterfactual joint effect"))) %>%
  group_by(sex) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  mutate(duration_percent = if_else(sim == "Observed data", "Reference",
                                    sprintf("-%1.1f%%", duration_percent))) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = duration, y =  reorder(sim, -duration)), 
           stat = "identity", fill = "dark gray") +
  geom_text(aes(x = duration, y =  reorder(sim, -duration),
                label = duration_percent),
            hjust = 1, size = 3) +
  facet_grid(sex~.) +
  labs(y = "Simulation*", x = "Educational difference in mean duration of MetS (years)") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0),
        axis.ticks.y = element_blank()); mean_duration_of_dif_plot_sens
ggsave("./plots/mean_duration_of_dif_plot_sens.pdf", plot = mean_duration_of_dif_plot_sens)


# merging various population-level parameters into one data frame and
# calculating percentage change of all parameters in counterfactuals
percent_diff_mets <- lc_prev_mets_both %>%
  relocate(lc_prev, .after = last_col()) %>%
  left_join(mean_onset_mets_both, by = c("sim", "sex", "analysis", "education")) %>%
  left_join(duration_of_mets_both, by = c("sim", "sex", "analysis", "education")) %>%
  group_by(sex, analysis) %>%
  mutate(across(where(is.numeric), 
                .fns = list(percent = ~ round((abs((.x - first(.))))/first(.) * 100, 1)))) %>%
  ungroup() %>%
  mutate(across(5:7, ~round(.,3))) %>%
  relocate(lc_prev_percent, .after = lc_prev) %>%
  relocate(mean_onset_percent, .after = mean_onset) %>%
  arrange(analysis, sex, -lc_prev_percent) 

write_xlsx(percent_diff_mets, path = "percent_diff_mets.xlsx")
