# load libraries ----------------------------------------------------------
library(MicSim)
library(tidyverse)
library(scales)

# setting path to data files
file_path <- "./data_files/" 

# setting path to store output
mc_path <- "./mc_files/"

# defining parameters used to repeat simulations with different N
max_iter <- 20
min_N <- 100000
max_N <- 500000 
step_size <- 100000
n_iter <- length(seq(from = min_N, to = max_N, by = step_size))
n_cores <- 24

# creating age variable values equal midpoint of our age groups
# for use in the parameterization
age <- c(19.5,seq(23.5,63.5,5))


# MetS parameterization ---------------------------------------------------

# loading observed data
df_mets <- read.csv(paste0(file_path, "rates_mets.csv")) %>%
  mutate(across(everything(), ~./100))

# estimating parameters to predict incidence rates 
model_mets <- list()
for (i in c(1:4)){
  y <- df_mets[,i]
  model_mets[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# estimating parameters to predict remittance rates 
model_mets.r <- list()
for (i in c(5:8)){
  y <- df_mets[,i]
  model_mets.r[[i]] <-  lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}


# creating functions that will be used to estimate transition
# probabilities

# incidence: males with low education 
noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(model_mets[[1]])[3]
  r <- coef(model_mets[[1]])[2]
  Po <- coef(model_mets[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# incidence: females with low education 
noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(model_mets[[2]])[3]
  r <- coef(model_mets[[2]])[2]
  Po <- coef(model_mets[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# incidence: males with high education
noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(model_mets[[3]])[3]
  r <- coef(model_mets[[3]])[2]
  Po <- coef(model_mets[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age))) 
  return(rate)
} 

# incidence: females with high education
noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(model_mets[[4]])[3]
  r <- coef(model_mets[[4]])[2]
  Po <- coef(model_mets[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age))) 
  return(rate)
} 


# remittance: males with low education
ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(model_mets.r[[5]])[[1]]
  a  <- coef(model_mets.r[[5]])[[2]]
  a2 <- coef(model_mets.r[[5]])[[3]]
  a3 <- coef(model_mets.r[[5]])[[4]]
  a4 <- coef(model_mets.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# remittance: females with low education
ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(model_mets.r[[6]])[[1]]
  a  <- coef(model_mets.r[[6]])[[2]]
  a2 <- coef(model_mets.r[[6]])[[3]]
  a3 <- coef(model_mets.r[[6]])[[4]]
  a4 <- coef(model_mets.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# remittance: males with high education 
ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(model_mets.r[[7]])[[1]]
  a  <- coef(model_mets.r[[7]])[[2]]
  a2 <- coef(model_mets.r[[7]])[[3]]
  a3 <- coef(model_mets.r[[7]])[[4]]
  a4 <- coef(model_mets.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# remittance: females with high education
ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(model_mets.r[[8]])[[1]]
  a  <- coef(model_mets.r[[8]])[[2]]
  a2 <- coef(model_mets.r[[8]])[[3]]
  a3 <- coef(model_mets.r[[8]])[[4]]
  a4 <- coef(model_mets.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# mortality rates ---------------------------------------------------------

#mortality rate is set to 0 in all simulations
mortRates <- function(age,calTime,duration){
  rate <- 0
  return(rate)
}


# MetS MC variance --------------------------------------------------------

# creating data frames to store results for life course prevalence, 
# mean age of onset, and mean duration of MetS

mc_prev_mets <- data.frame(matrix(ncol = n_iter + 3, nrow = max_iter*4))
colnames(mc_prev_mets) <- c(as.character(1:n_iter), "iteration", "education", "gender") 
mc_prev_mets$iteration <- rep(1:max_iter, each = 4)
mc_prev_mets$education <- rep(c("High", "High", "Low", "Low"), max_iter)
mc_prev_mets$gender <- rep(c("Female", "Male", "Female", "Male"), max_iter)

mc_age_mets <- data.frame(matrix(ncol = n_iter + 3, nrow = max_iter*4))
colnames(mc_age_mets) <- c(as.character(1:n_iter), "iteration", "education", "gender") 
mc_age_mets$iteration <- rep(1:max_iter, each = 4)
mc_age_mets$education <- rep(c("High", "High", "Low", "Low"), max_iter)
mc_age_mets$gender <- rep(c("Female", "Male", "Female", "Male"), max_iter)

mc_time_mets <- data.frame(matrix(ncol = n_iter + 3, nrow = max_iter*4))
colnames(mc_time_mets) <- c(as.character(1:n_iter), "iteration", "education", "gender") 
mc_time_mets$iteration <- rep(1:max_iter, each = 4)
mc_time_mets$education <- rep(c("High", "High", "Low", "Low"), max_iter)
mc_time_mets$gender <- rep(c("Female", "Male", "Female", "Male"), max_iter)

# setting variable to iterate through columns of df that store results
# each column represents a different sample size
col_num <- 1

# for loop to iterate through desired sample sizes
for (N in seq(from = min_N, to = max_N, by = step_size)) {
  
  # defining initial pop 
  # birth cohort 2000 (they enter simulation in 2018 at age 18)
  set.seed(20210607)
  initBirthDatesRange <- chron(dates=c("1/1/2000","31/12/2000"), format=c(dates="d/m/Y"),out.format=c(dates="d/m/year"))
  birthDates <- dates(initBirthDatesRange[1] + runif(N, min=0, max=diff(initBirthDatesRange)))
  
  # setting minimum and maximum age 
  minage <- 18
  maxage <- 65
  
  # simulation horizon
  simHorizon <- setSimHorizon(startDate="01/01/2018", endDate="31/12/2065")
  
  # initial states: intersection between gender, education level and health status
  states_mets <- c("noms.lo.m","noms.lo.f","noms.hi.m","noms.hi.f",
                   "ms.lo.m","ms.lo.f","ms.hi.m","ms.hi.f")
  
  # absorbing states: dead (must be included even though mortality rate set to zero)
  absStates <- "dead"
  
  # all states 
  stateSpace_mets <- expand.grid(states_mets)
  
  # prevalence of MetS in starting pop depends on a stochastic process driven by 
  # percentages provided by regressions. assume the sample has an equal distribution 
  # of males and females with low and high education  
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
  
  # setting variable to iterate through rows of df storing results
  # starts off as 1:4 because there are 4 groups per iteration:
  # females with high education, males with high education, 
  # females with low education, and males with low education
  mc_rows <- 1:4
  
  # for loop to repeat simulation with a set N. 
  # max_iter determines the number of times the analysis is run per N
  for (iter in 1:max_iter) {
    # running simulations -----------------------------------------------------
    seed = N + 100 * iter
    simpop_mets <- micSimParallel(initPop=initPop_mets, transitionMatrix=transitionMatrix_mets, absStates=absStates,
                                  maxAge=maxage,simHorizon=simHorizon, cores = n_cores, seeds = seed) 
    
    # converting output to spell format 
    # first duplicating data
    msm_mets <- simpop_mets %>% arrange(ID, transitionTime)
    
    # creating index to identify each transition
    msm_mets <- msm_mets %>% group_by(ID) %>% mutate(index=1:n()) %>% ungroup()
    
    # only selecting first row per individual and adding labels for 
    # gender, education and health status at start of model
    obs1_mets <- msm_mets %>%
      filter(index == 1) %>%
      mutate(gender = if_else(str_detect(initState, ".m$"), "Male", "Female"), 
             education = if_else(str_detect(initState, ".lo."), "Low", "High"),
             health = if_else(str_detect(initState, "noms."), "NoMS", "MS"),
             agestart = minage)
    
    # recoding health status 
    msm_mets$health <-  if_else(str_detect(msm_mets$To, "noms."), "NoMS", "MS")
    
    # rearrange data from msm format to spell format (one line per state)
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
    
    # assign index 
    spell_mets <- spell_mets %>% group_by(ID) %>% mutate(index=1:n()) %>% ungroup()
    
    # estimating the required metrics 
    mc_prev_mets[mc_rows,col_num] <- spell_mets %>%
      mutate(health_num = if_else(health == "MS", 1, 0)) %>%
      # arranges so that instances of MetS appear first
      arrange(health) %>%
      # only keeps first row per individual
      distinct(ID, .keep_all = T) %>%
      group_by(education, gender) %>%
      summarise(prev = mean(health_num) * 100, .groups = "drop")  %>%
      select(prev)
    
    mc_age_mets[mc_rows,col_num] <- spell_mets %>%
      # if no MetS when entering cohort then first instance
      # would be recorded at index 2
      filter(index == 2, health == "MS") %>%
      group_by(education, gender) %>%
      summarise(age_onset = mean(agestart), .groups = "drop") %>%
      select(age_onset)
    
    mc_time_mets[mc_rows,col_num] <- spell_mets %>%
      filter(health == "MS") %>%
      group_by(education, gender) %>%
      summarise(time_with = sum(length)/length(unique(ID)), .groups = "drop") %>% 
      select(time_with)
    
    mc_rows <- mc_rows + 4
  }
  
  col_num <- col_num + 1
}

# exporting df, so I can analyze results later
save(mc_prev_mets, file = paste0(mc_path, "mc_prev_mets.rda"))
save(mc_age_mets, file = paste0(mc_path, "mc_age_mets.rda"))
save(mc_time_mets, file = paste0(mc_path, "mc_time_mets.rda"))

# Analyzing the monte carlo variance --------------------------------------

# loading the files
load(paste0(mc_path, "mc_prev_mets.rda"))
load(paste0(mc_path, "mc_age_mets.rda"))
load(paste0(mc_path, "mc_time_mets.rda"))

# calculating the standard deviation of the population-level parameters
# across the various sample sizes, and then plotting this
mc_prev_mets %>%
  group_by(education, gender) %>%
  summarise(across(1:n_iter, ~sd(.)), .groups = "drop") %>%
  pivot_longer(3:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x = as.numeric(name) * step_size, y = value, group = education, color = education)) +
  facet_grid(gender ~ .) +
  scale_x_continuous(labels = comma) +
  labs(y = "Standard deviation*", x = "Sample size", 
       colour = "Education")
ggsave(filename = "MetS MC variance prev.pdf", path = "./plots", 
       plot = last_plot(), width = 7, height = 5, unit = "in")

mc_age_mets %>%
  group_by(education, gender) %>%
  summarise(across(1:n_iter, ~sd(.))) %>%
  pivot_longer(3:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x = as.numeric(name) * step_size, y = value, group = education, color = education)) +
  facet_grid(gender ~ .) +
  scale_x_continuous(labels = comma) +
  labs(y = "Standard deviation*", x = "Sample size", 
       #title = "MetS: Monte carlo variance for mean age of onset",
       colour = "Education")
ggsave(filename = "MetS MC variance onset.pdf", path = "./plots", 
       plot = last_plot(), width = 7, height = 5, unit = "in")

mc_time_mets %>%
  group_by(education, gender) %>%
  summarise(across(1:n_iter, ~sd(.))) %>%
  pivot_longer(3:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x = as.numeric(name) * step_size, y = value, group = education, color = education)) +
  facet_grid(gender ~ .) +
  scale_x_continuous(labels = comma) +
  labs(y = "Standard deviation*", x = "Sample size", 
       #title = "MetS: Monte carlo variance for meant time spent",
       colour = "Education")
ggsave(filename = "MetS MC variance duration.pdf", path = "./plots", 
       plot = last_plot(), width = 7, height = 5, unit = "in")



