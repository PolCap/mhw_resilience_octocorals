# ---------------------------------------------------------------------------- #
# - FILE NAME:   DemoFunctions.R         
# - DATE:        14/07/2023
# - DESCRIPTION: Functions for the manuscript Marine heatwaves select for resistant 
#               Mediterranean octocoral populations at the expense speed of recovery
# - AUTHORS:     Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #

# These fucntions are required to produce some of the results. To run these 
# functions at least the packages matrixcalc, popdemo, and popbio need to be 
# installed. Also, tidyverse is used in some parts but is not essential. 

# Vital rates function ---------------------------------------------------------

# This is a generic function created to estimate different vital rates from 
# demographic data collected from permanent transects, with size measurements. 
# It is only applicable to populations where individuals were tagged and measured
# in size in yearly transitions (a yearly monitoring). 
# The function does different things: 
# - Classifies the different sizes into differnt classes that will be 
# pre-defined by the author.
# - Estimates the survival of organisms by accounting for the number of
# individuals alive at time t+1 divided by the total of individuals at time t.
# - Estimates the transition of individuals from time t to time t+1 into
# the next size class. 
# - Calculates the shrinkage of individuals into lower size classes at
# time t+1. Here the function calculates considering up to a drop into three
# size classes.  
# - Calculates the reproductive output based on the fertility of the 
# individuals and the number of recruits per adult. 

vital_rates <- function(size, size_next, size_next2,
                        class, breaks, fert, as_list=TRUE,
                        ad_first=TRUE){
  
  # Create a tibble with the sizes
  
  if(sum(size_next2, na.rm=T)==0){
    data <- tibble(size=size, size_next=size_next)
  }else{
    data <- tibble(size=size, size_next=size_next,size_next2=size_next2)
  }
  
  # Create a tibble with the fertility (ordered by the size classes)
  
  f <- tibble(class=class, fert=fert)
  
  # Classify into size classes
  
  data <- data %>% 
    mutate(across(starts_with("size"),
                  ~cut(.x, 
                       breaks = breaks,
                       labels = class,
                       ordered_result = TRUE,
                       right= TRUE),
                  .names = "class_{.col}")) %>% 
    # Remove those individuals that are not recruits and where not considered in
    # the previous year
    filter(
      if ("class_size_next2" %in% names(.)) {
        !is.na(class_size) | class_size_next2 == class[2]
      } else {
        !is.na(class_size)
      })
  
  # Count number of recruits in pos year
  
  suppressWarnings(num_recr <- sum(is.na(data$class_size_next) & data$class_size_next2 == class[2], 
                                   na.rm = TRUE))
  
  if(sum(size_next2, na.rm = T)!=0){
    # Calculate survival of the first class 
    
    surv_first <- data %>%
      filter(!is.na(class_size_next)) %>%
      group_by(class_size_next) %>%
      summarise(T = n(),
                S = sum(!is.na(class_size_next2))) %>%
      mutate(s = S / T, 0) %>% 
      filter(class_size_next==class[2]) %>% 
      select(s)
    
    # Estimate the past recruits
    
    num_recr <- num_recr/as.numeric(surv_first$s)  
  } else{
    num_recr
  }
  
  # Calculate survival, growth, and shrinkage
  data <- data %>% 
    filter(!is.na(class_size)) %>% 
    mutate(# Survival
      survival=ifelse(is.na(size_next)|size_next==0, 0,1),
      # Growth 
      growth=ifelse(class_size<class_size_next, 1, 0),
      # Shrinkage
      shrinkage_one=ifelse(as.numeric(class_size)-as.numeric(class_size_next)==1, 
                           1, 0),
      shrinkage_two=ifelse(as.numeric(class_size)-as.numeric(class_size_next)==2, 
                           1, 0),
      shrinkage_three=ifelse(as.numeric(class_size)-as.numeric(class_size_next)==3, 
                             1, 0)) %>% 
    # Reproduction 
    # Define fertility by classes 
    left_join(f, by = c("class_size"="class")) %>% 
    mutate(fertility = size^2.61*fert) %>% 
    # Add all the values  
    group_by(class_size) %>% 
    summarise(T = n(),
              S = sum(survival, na.rm = T),
              G = sum(growth, na.rm = TRUE),  
              H = sum(shrinkage_one, na.rm = TRUE),
              HH = sum(shrinkage_two, na.rm = TRUE),
              HHH = sum(shrinkage_three, na.rm = TRUE), 
              fec=sum(fertility, na.rm = TRUE),
              num_recr=num_recr) %>% 
    mutate(s = S / T, 0,
           g = G / S,
           h = (H + HH)/(S - G),
           hl = (HH)/(H + HH)) %>% 
    ungroup() %>% 
    mutate(f=ifelse(num_recr > 0, num_recr * (fec/sum(fec))/T, 0)) %>% 
    # if Na or NaN, change to 0:
    mutate(across(-class_size,
                  ~ifelse(is.na(.x) | is.nan(.x), 0, .x)),
           #Change class_size into an ordered factor
           class_size=ordered(class_size, class)) %>% 
    select(class_size, f, s, g, h, hl) %>% 
    gather(key = "vital_rate" , 
           value = "value",f,s,g,h,hl) %>% 
    mutate(class_size=fct_recode(class_size,
                                 "1" = "0-0.3",
                                 "2" = "0.3-3", 
                                 "3" = "3-10", 
                                 "4" = "11-20",
                                 "5" = "21-30",
                                 "6" = "31-40",
                                 "7" = ">40"),
           vital_rate=paste0(vital_rate, class_size)) %>% 
    select(-class_size)
  
  if(ad_first==TRUE & sum(data$value)>0){
    data <- data %>% rbind(data.frame(vital_rate=c("s1","g1"), 
                                      value=c(data$value[data$vital_rate=="s2"], 1)))
  }else{
    data 
  }
  
  if(as_list==TRUE){
    
    data <-  data %>% split(f = .$vital_rate) %>% map("value")
    
  }
  
  return(data) 
  
}

# Function to decompose the matrix into demographic processes ------------------

# First, a separate function to create matF and matU

mat_decomp <- function(x){
  matF <- matU <- x
  matF[-1,] <- 0
  matU[1,] <- 0
  return(tibble(matF, matU))
}

# Second, function that returns the different components

mat_comp <- function(mat) {
  # decompose the matrix into matF and matU
  mats<- mat_decomp(mat)
  # Contributions to the reproduction
  fecundity <- sum(mats$matF)
  # Contributions to stasis
  stasis <- sum(diag(mats$matU))
  # Contributions to shrinkage
  shrinkage <- sum(mats$matU[upper.tri(mats$matU)])
  # Contributions to growth
  growth <- sum(mats$matU[lower.tri(mats$matU)])
  
  return(data.frame(fecundity, growth,
                    shrinkage, stasis))
}

# Functions to estimate resilience components elasticities --------------------- 
# (code modified from Cant et al. (2023) EcoLetts) 

# Damping ratio sensitivity and elasticity

elas_dr <- function(mat, components=T) {
  # Calculate damping ratio
  dr <- popdemo::dr(mat)
  # Calculate sensitivity for dominant eigenvalue
  s1 <- sens(mat, eval = "max", all = FALSE)
  
  # Calculate sensitivity for subdominant eigenvalue
  s2 <- sens(mat, eval = 2, all = FALSE)
  
  lambda2 <- abs(eigen(mat)$values)[2]
  
  # Calculate element-level sensitivities of the damping ratio
  srho <- (1 / lambda2) * (s1 - (dr / lambda2) * (Re(s2) + Im(s2)))
  
  # Calculate damping ratio elasticity
  erho <- matrixcalc::hadamard.prod((1/dr)*srho, mat)

    # Decompose the matrix into its components
  
  comp <- mat_comp(erho)
  
  # Return depending on the specifics of the function
  
  if(components==T){
    return(comp)
  }else{
    return(tibble(srho, erho))
  }
}

# Sensitivity and Elasticity of first time step attenuation using the brute force 
# method by Morris and Doak (2002)

elas_res <- function(mat, perturb = 0.01,
                       components=T) {
  
  # Calculate initial first step attenuation
  atten1 <- 1-popdemo::reac(mat, bound = "lower")
  
  drow <- dim(mat)[1]
  dcol <- dim(mat)[2]
  sens_a <- matrix(NA, drow, dcol)
  elas_a <- matrix(NA, drow, dcol)
  
  for (i in 1:drow) {
    for (j in 1:dcol) {
      mat_new <- mat
      mat_new[i, j] <- mat_new[i, j] + (mat_new[i, j] * perturb)
      element_diff <- mat_new[i, j] * perturb
      
      # Calculatw first step attenuation after perturbation
      atten2 <- popdemo::reac(mat_new, bound = "lower")
      
      # Calculate element-level sensitivity
      sens_a[i, j] <- (atten2 - atten1) / element_diff
      
      # Calculate element-level elasticity
      elas_a[i, j] <- ((atten2 - atten1) / atten1) / (element_diff / mat[i, j])
    }
  }
  
  # Replace NaN and infinite values with 0
  sens_a[is.na(sens_a) | is.infinite(sens_a)] <- 0
  elas_a[is.na(elas_a) | is.infinite(elas_a)] <- 0
  
  # Decompose the elasticity matrix
  comp <- mat_comp(elas_a)
  
  if(isTRUE(components)){
    return(comp)
  }else{
    return(tibble(sens_a, elas_a))
  }
}

elas_rec <- function(mat, perturb = 0.01,
                     components=T) {
  
  # Calculate initial first step attenuation
  atten1 <- popdemo::dr(mat)
  
  drow <- dim(mat)[1]
  dcol <- dim(mat)[2]
  sens_a <- matrix(NA, drow, dcol)
  elas_a <- matrix(NA, drow, dcol)
  
  for (i in 1:drow) {
    for (j in 1:dcol) {
      mat_new <- mat
      mat_new[i, j] <- mat_new[i, j] + (mat_new[i, j] * perturb)
      element_diff <- mat_new[i, j] * perturb
      
      # Calculatw first step attenuation after perturbation
      atten2 <- popdemo::dr(mat_new)
      
      # Calculate element-level sensitivity
      sens_a[i, j] <- (atten2 - atten1) / element_diff
      
      # Calculate element-level elasticity
      elas_a[i, j] <- ((atten2 - atten1) / atten1) / (element_diff / mat[i, j])
    }
  }
  
  # Replace NaN and infinite values with 0
  sens_a[is.na(sens_a) | is.infinite(sens_a)] <- 0
  elas_a[is.na(elas_a) | is.infinite(elas_a)] <- 0
  
  # Decompose the elasticity matrix
  comp <- mat_comp(elas_a)
  
  if(isTRUE(components)){
    return(comp)
  }else{
    return(tibble(sens_a, elas_a))
  }
}

# Function to project the population estimating population quasi-extinction, 
# resilience metrics, and elasticies to 

pop_proj <- function(mats, mats_mhw, prob, n0, tmax, reps) {
  
  # Generate empty data frame
  
  data <-  expand_grid(rep=seq(1, reps, by=1),
                       time=seq(1, tmax, by=1),
                       prob=prob,
                       resistance = NA, 
                       recovery = NA,
                       structure=NA, 
                       n=NA)
  data <- as_tibble(data)
  data$structure[data$time==1] <- list(n0)
  data$n[data$time==1] <- sum(n0)
  
  for(r in 1:reps){
    # Create a for loop to project over time 
    for (t in 2:tmax) {
      u <- rbinom(1, 1, prob) # Binary to determine whether there is MHW
      A <- ifelse(u == 1, sample(mats_mhw, 1), sample(mats, 1)) # Sample the matrices
      
      # Calculate the resilience components
      
      data$resistance[data$time==t-1& data$rep==r] <- reac(A[[1]], bound = "lower")
      data$recovery[data$time==t-1& data$rep==r] <- dr(A[[1]])

      # Estimate the abundance
      data$structure[data$time==t & data$rep==r] <- list(A[[1]] %*% data$structure[data$time==t-1 & data$rep==r][[1]])
      data$n[data$time==t & data$rep==r] <- sum(A[[1]] %*% data$structure[data$time==t-1 & data$rep==r][[1]])
      }
    }
  return(data)
}
