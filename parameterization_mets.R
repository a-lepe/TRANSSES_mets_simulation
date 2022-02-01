# load libraries and read in data -----------------------------------------
library(tidyverse) 
#library(measures)

# setting path to data files
file_path <- "./data_files/" 

# creating age variable values equal midpoint of our age groups
age <- c(19.5,seq(23.5,63.5,5))

# observed data
df_mets <- read.csv(paste0(file_path, "rates_mets.csv")) %>%
  mutate(across(everything(), ~./100))

# counterfactual data
cf_df_mets_smk <- read.csv(paste0(file_path, "rates_mets_smk.csv")) %>%
  mutate(across(everything(), ~./100))

cf_df_mets_alc <- read.csv(paste0(file_path, "rates_mets_alc.csv")) %>%
  mutate(across(everything(), ~./100))

cf_df_mets_diet <- read.csv(paste0(file_path, "rates_mets_diet.csv")) %>%
  mutate(across(everything(), ~./100))

cf_df_mets_hl <- read.csv(paste0(file_path, "rates_mets_hl.csv")) %>%
  mutate(across(everything(), ~./100))

cf_df_mets_joint <- read.csv(paste0(file_path, "rates_mets_joint.csv")) %>%
  mutate(across(everything(), ~./100))


# estimate parameters to predict transition rates -------------------------
# observed data mets: from not having to having MS: logistic
model_mets <- list()
for (i in c(1:4)){
  y <- df_mets[,i]
  model_mets[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# observed data mets: from having to not having MS: quadratic
model_mets.r <- list()
for (i in c(5:8)){
  y <- df_mets[,i]
  model_mets.r[[i]] <-  lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}

# counterfactual data mets smk: from not having to having MS: logistic
cf_model_smk <- list()
for (i in c(1:4)){
  y <- cf_df_mets_smk[,i]
  cf_model_smk[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# counterfactual data mets smk: having to not having MS: quadratic 
cf_model_smk.r <- list()
for (i in c(5:8)){
  y <- cf_df_mets_smk[,i]
  cf_model_smk.r[[i]] <- lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}

# counterfactual data mets alc: from not having to having MS: logistic 
cf_model_alc <- list()
for (i in c(1:4)){
  y <- cf_df_mets_alc[,i]
  cf_model_alc[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# counterfactual data mets alc: having to not having MS: quadratic
cf_model_alc.r <- list()
for (i in c(5:8)){
  y <- cf_df_mets_alc[,i]
  cf_model_alc.r[[i]] <- lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}

# counterfactual data mets diet: from not having to having MS: logistic
cf_model_diet <- list()
for (i in c(1:4)){
  y <- cf_df_mets_diet[,i]
  cf_model_diet[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# counterfactual data mets diet: having to not having MS: quadratic
cf_model_diet.r <- list()
for (i in c(5:8)){
  y <- cf_df_mets_diet[,i]
  cf_model_diet.r[[i]] <- lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}

# counterfactual data mets hl: from not having to having MS: logistic
cf_model_hl <- list()
for (i in c(1:4)){
  y <- cf_df_mets_hl[,i]
  cf_model_hl[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# counterfactual data mets hl: having to not having MS: quadratic
cf_model_hl.r <- list()
for (i in c(5:8)){
  y <- cf_df_mets_hl[,i]
  cf_model_hl.r[[i]] <- lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}

# counterfactual data mets joint: from not having to having MS: logistic
cf_model_joint <- list()
for (i in c(1:4)){
  y <- cf_df_mets_joint[,i]
  cf_model_joint[[i]] <- nls(y ~ K / (1 + exp(Po + r *age)), start=list(Po=3, r=-.1, K=0.19))
}

# counterfactual data mets joint: having to not having MS: quadratic
cf_model_joint.r <- list()
for (i in c(5:8)){
  y <- cf_df_mets_joint[,i]
  cf_model_joint.r[[i]] <- lm(y ~ age + I(age^2) + I(age^3) + I(age^4))
}


# functions to estimate transition rates using model coefficients ---------
# calTime and duration are needed even though we don't use those parameters

# functions using observed data MetS  -------------------------------------
# observed data MetS: from not having MS to having MS
# males with low education
noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(model_mets[[1]])[3]
  r <- coef(model_mets[[1]])[2]
  Po <- coef(model_mets[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with low education
noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(model_mets[[2]])[3]
  r <- coef(model_mets[[2]])[2]
  Po <- coef(model_mets[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# males with high education
noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(model_mets[[3]])[3]
  r <- coef(model_mets[[3]])[2]
  Po <- coef(model_mets[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age))) 
  return(rate)
} 

# females with high education
noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(model_mets[[4]])[3]
  r <- coef(model_mets[[4]])[2]
  Po <- coef(model_mets[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age))) 
  return(rate)
} 

# oberved data MetS: from having MS to not having MS
# males with low education
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

# females with low education
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

# males with high education
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

# females with high education
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


# functions using counterfactual data MetS smk  ---------------------------
# counterfactual data MetS smk: from not having MS to having MS
# males with low education
cf_smk_noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_smk[[1]])[3]
  r <- coef(cf_model_smk[[1]])[2]
  Po <- coef(cf_model_smk[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# females with low education
cf_smk_noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_smk[[2]])[3]
  r <- coef(cf_model_smk[[2]])[2]
  Po <- coef(cf_model_smk[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# males with high education
cf_smk_noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_smk[[3]])[3]
  r <- coef(cf_model_smk[[3]])[2]
  Po <- coef(cf_model_smk[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with high education
cf_smk_noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_smk[[4]])[3]
  r <- coef(cf_model_smk[[4]])[2]
  Po <- coef(cf_model_smk[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))    
  return(rate)
} 

# counterfactual data MetS smk: from having MS to not having MS
# males with low education
cf_smk_ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(cf_model_smk.r[[5]])[[1]]
  a  <- coef(cf_model_smk.r[[5]])[[2]]
  a2 <- coef(cf_model_smk.r[[5]])[[3]]
  a3 <- coef(cf_model_smk.r[[5]])[[4]]
  a4 <- coef(cf_model_smk.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with low education
cf_smk_ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(cf_model_smk.r[[6]])[[1]]
  a  <- coef(cf_model_smk.r[[6]])[[2]]
  a2 <- coef(cf_model_smk.r[[6]])[[3]]
  a3 <- coef(cf_model_smk.r[[6]])[[4]]
  a4 <- coef(cf_model_smk.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# males with high education
cf_smk_ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(cf_model_smk.r[[7]])[[1]]
  a  <- coef(cf_model_smk.r[[7]])[[2]]
  a2 <- coef(cf_model_smk.r[[7]])[[3]]
  a3 <- coef(cf_model_smk.r[[7]])[[4]]
  a4 <- coef(cf_model_smk.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with high education
cf_smk_ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(cf_model_smk.r[[8]])[[1]]
  a  <- coef(cf_model_smk.r[[8]])[[2]]
  a2 <- coef(cf_model_smk.r[[8]])[[3]]
  a3 <- coef(cf_model_smk.r[[8]])[[4]]
  a4 <- coef(cf_model_smk.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# functions using counterfactual data MetS alc ----------------------------
# counterfactual data MetS alc: from not having MS to having MS
# males with low education
cf_alc_noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_alc[[1]])[3]
  r <- coef(cf_model_alc[[1]])[2]
  Po <- coef(cf_model_alc[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# females with low education
cf_alc_noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_alc[[2]])[3]
  r <- coef(cf_model_alc[[2]])[2]
  Po <- coef(cf_model_alc[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# males with high education
cf_alc_noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_alc[[3]])[3]
  r <- coef(cf_model_alc[[3]])[2]
  Po <- coef(cf_model_alc[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with high education
cf_alc_noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_alc[[4]])[3]
  r <- coef(cf_model_alc[[4]])[2]
  Po <- coef(cf_model_alc[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))    
  return(rate)
} 

# counterfactual data MetS alc: from having MS to not having MS
# males with low education
cf_alc_ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(cf_model_alc.r[[5]])[[1]]
  a  <- coef(cf_model_alc.r[[5]])[[2]]
  a2 <- coef(cf_model_alc.r[[5]])[[3]]
  a3 <- coef(cf_model_alc.r[[5]])[[4]]
  a4 <- coef(cf_model_alc.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with low education
cf_alc_ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(cf_model_alc.r[[6]])[[1]]
  a  <- coef(cf_model_alc.r[[6]])[[2]]
  a2 <- coef(cf_model_alc.r[[6]])[[3]]
  a3 <- coef(cf_model_alc.r[[6]])[[4]]
  a4 <- coef(cf_model_alc.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# males with high education
cf_alc_ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(cf_model_alc.r[[7]])[[1]]
  a  <- coef(cf_model_alc.r[[7]])[[2]]
  a2 <- coef(cf_model_alc.r[[7]])[[3]]
  a3 <- coef(cf_model_alc.r[[7]])[[4]]
  a4 <- coef(cf_model_alc.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with high education
cf_alc_ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(cf_model_alc.r[[8]])[[1]]
  a  <- coef(cf_model_alc.r[[8]])[[2]]
  a2 <- coef(cf_model_alc.r[[8]])[[3]]
  a3 <- coef(cf_model_alc.r[[8]])[[4]]
  a4 <- coef(cf_model_alc.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# functions using counterfactual data MetS diet ---------------------------
# counterfactual data MetS diet: from not having MS to having MS
# males with low education
cf_diet_noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_diet[[1]])[3]
  r <- coef(cf_model_diet[[1]])[2]
  Po <- coef(cf_model_diet[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# females with low education
cf_diet_noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_diet[[2]])[3]
  r <- coef(cf_model_diet[[2]])[2]
  Po <- coef(cf_model_diet[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# males with high education
cf_diet_noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_diet[[3]])[3]
  r <- coef(cf_model_diet[[3]])[2]
  Po <- coef(cf_model_diet[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with high education
cf_diet_noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_diet[[4]])[3]
  r <- coef(cf_model_diet[[4]])[2]
  Po <- coef(cf_model_diet[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))    
  return(rate)
} 

# counterfactual data MetS diet: from having MS to not having MS
# males with low education
cf_diet_ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(cf_model_diet.r[[5]])[[1]]
  a  <- coef(cf_model_diet.r[[5]])[[2]]
  a2 <- coef(cf_model_diet.r[[5]])[[3]]
  a3 <- coef(cf_model_diet.r[[5]])[[4]]
  a4 <- coef(cf_model_diet.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with low education
cf_diet_ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(cf_model_diet.r[[6]])[[1]]
  a  <- coef(cf_model_diet.r[[6]])[[2]]
  a2 <- coef(cf_model_diet.r[[6]])[[3]]
  a3 <- coef(cf_model_diet.r[[6]])[[4]]
  a4 <- coef(cf_model_diet.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# males with high education
cf_diet_ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(cf_model_diet.r[[7]])[[1]]
  a  <- coef(cf_model_diet.r[[7]])[[2]]
  a2 <- coef(cf_model_diet.r[[7]])[[3]]
  a3 <- coef(cf_model_diet.r[[7]])[[4]]
  a4 <- coef(cf_model_diet.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with high education
cf_diet_ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(cf_model_diet.r[[8]])[[1]]
  a  <- coef(cf_model_diet.r[[8]])[[2]]
  a2 <- coef(cf_model_diet.r[[8]])[[3]]
  a3 <- coef(cf_model_diet.r[[8]])[[4]]
  a4 <- coef(cf_model_diet.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# functions using counterfactual data MetS hl -----------------------------
# counterfactual data MetS hl: from not having MS to having MS
# males with low education
cf_hl_noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_hl[[1]])[3]
  r <- coef(cf_model_hl[[1]])[2]
  Po <- coef(cf_model_hl[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# females with low education
cf_hl_noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_hl[[2]])[3]
  r <- coef(cf_model_hl[[2]])[2]
  Po <- coef(cf_model_hl[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# males with high education
cf_hl_noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_hl[[3]])[3]
  r <- coef(cf_model_hl[[3]])[2]
  Po <- coef(cf_model_hl[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with high education
cf_hl_noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_hl[[4]])[3]
  r <- coef(cf_model_hl[[4]])[2]
  Po <- coef(cf_model_hl[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))    
  return(rate)
} 

# counterfactual data MetS hl: from having MS to not having MS
# males with low education
cf_hl_ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(cf_model_hl.r[[5]])[[1]]
  a  <- coef(cf_model_hl.r[[5]])[[2]]
  a2 <- coef(cf_model_hl.r[[5]])[[3]]
  a3 <- coef(cf_model_hl.r[[5]])[[4]]
  a4 <- coef(cf_model_hl.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with low education
cf_hl_ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(cf_model_hl.r[[6]])[[1]]
  a  <- coef(cf_model_hl.r[[6]])[[2]]
  a2 <- coef(cf_model_hl.r[[6]])[[3]]
  a3 <- coef(cf_model_hl.r[[6]])[[4]]
  a4 <- coef(cf_model_hl.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# males with high education
cf_hl_ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(cf_model_hl.r[[7]])[[1]]
  a  <- coef(cf_model_hl.r[[7]])[[2]]
  a2 <- coef(cf_model_hl.r[[7]])[[3]]
  a3 <- coef(cf_model_hl.r[[7]])[[4]]
  a4 <- coef(cf_model_hl.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with high education
cf_hl_ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(cf_model_hl.r[[8]])[[1]]
  a  <- coef(cf_model_hl.r[[8]])[[2]]
  a2 <- coef(cf_model_hl.r[[8]])[[3]]
  a3 <- coef(cf_model_hl.r[[8]])[[4]]
  a4 <- coef(cf_model_hl.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# functions using counterfactual data MetS joint --------------------------
# counterfactual data MetS joint: from not having MS to having MS
# males with low education
cf_joint_noms_ms.lo.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_joint[[1]])[3]
  r <- coef(cf_model_joint[[1]])[2]
  Po <- coef(cf_model_joint[[1]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# females with low education
cf_joint_noms_ms.lo.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_joint[[2]])[3]
  r <- coef(cf_model_joint[[2]])[2]
  Po <- coef(cf_model_joint[[2]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))   
  return(rate)
} 

# males with high education
cf_joint_noms_ms.hi.m <- function(age,calTime,duration){
  
  K <- coef(cf_model_joint[[3]])[3]
  r <- coef(cf_model_joint[[3]])[2]
  Po <- coef(cf_model_joint[[3]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))  
  return(rate)
} 

# females with high education
cf_joint_noms_ms.hi.f <- function(age,calTime,duration){
  
  K <- coef(cf_model_joint[[4]])[3]
  r <- coef(cf_model_joint[[4]])[2]
  Po <- coef(cf_model_joint[[4]])[1]
  
  rate <- ifelse(age < 18, 0, K / (1 + exp(Po + r *age)))    
  return(rate)
} 

# counterfactual data MetS joint: from having MS to not having MS
# males with low education
cf_joint_ms_noms.lo.m <- function(age,calTime,duration){
  b  <- coef(cf_model_joint.r[[5]])[[1]]
  a  <- coef(cf_model_joint.r[[5]])[[2]]
  a2 <- coef(cf_model_joint.r[[5]])[[3]]
  a3 <- coef(cf_model_joint.r[[5]])[[4]]
  a4 <- coef(cf_model_joint.r[[5]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with low education
cf_joint_ms_noms.lo.f <- function(age,calTime,duration){
  b  <- coef(cf_model_joint.r[[6]])[[1]]
  a  <- coef(cf_model_joint.r[[6]])[[2]]
  a2 <- coef(cf_model_joint.r[[6]])[[3]]
  a3 <- coef(cf_model_joint.r[[6]])[[4]]
  a4 <- coef(cf_model_joint.r[[6]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# males with high education
cf_joint_ms_noms.hi.m <- function(age,calTime,duration){
  b  <- coef(cf_model_joint.r[[7]])[[1]]
  a  <- coef(cf_model_joint.r[[7]])[[2]]
  a2 <- coef(cf_model_joint.r[[7]])[[3]]
  a3 <- coef(cf_model_joint.r[[7]])[[4]]
  a4 <- coef(cf_model_joint.r[[7]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 

# females with high education
cf_joint_ms_noms.hi.f <- function(age,calTime,duration){
  b  <- coef(cf_model_joint.r[[8]])[[1]]
  a  <- coef(cf_model_joint.r[[8]])[[2]]
  a2 <- coef(cf_model_joint.r[[8]])[[3]]
  a3 <- coef(cf_model_joint.r[[8]])[[4]]
  a4 <- coef(cf_model_joint.r[[8]])[[5]]
  
  rate <- ifelse(age < 18, 0, 
                 b + a*age + a2*age^2 + a3*age^3 + a4*age^4)
  return(rate)
} 


# mortality rates ---------------------------------------------------------
# mortality rate is set to 0
mortRates <- function(age,calTime,duration){
  rate <- 0
  return(rate)
}


# checking model fit ------------------------------------------------------
# observed data MetS
observed_mets <- c(df_mets[,1], df_mets[,2], df_mets[,3],
                   df_mets[,4], df_mets[,5], df_mets[,6],
                   df_mets[,7], df_mets[,8])
fit_mets <- c(noms_ms.lo.m(age), noms_ms.lo.f(age),
              noms_ms.hi.m(age), noms_ms.hi.f(age),
              ms_noms.lo.m(age), ms_noms.lo.f(age),
              ms_noms.hi.m(age), ms_noms.hi.f(age))
sex <- c("Male","Female")
education <- c("Low education","Low education","High education","High education")
transition_mets <- c(rep("Incidence",4),rep("Remittance",4))
check_mets <- data.frame(age=age, observed=observed_mets, fit=fit_mets,
                         sex=rep(sex,each=10), 
                         education=rep(education,each=10),
                         transition=rep(transition_mets,each=10)) %>%
  group_by(sex, education, transition) %>%
  mutate(mae = mean(abs(observed-fit))) %>%
  ungroup()

# calculating the mean absolute error
check_mets %>%
  group_by(transition, sex, education) %>%
  summarise(mae = mean(abs(observed-fit))) %>%
  ungroup()

# plotting observed vs expected values
ggplot(check_mets)+
  geom_point(aes(age,observed,group=sex,color=sex))+
  geom_line(aes(age,fit,group=sex,color=sex))+
  ylab("transition rate") +
  xlab("Age (years)") +
  facet_grid(education~transition)
ggsave("./plots/publications/transition_fit_mets.pdf", plot = last_plot())


rm(observed_mets, fit_mets, check_mets, transition_mets, 
   sex, education)
