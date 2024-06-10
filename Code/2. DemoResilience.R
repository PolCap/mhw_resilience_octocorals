# --------------------------------------------------------------------------------------- #
# - FILE NAME:   DemoResilience.R         
# - DATE:        20/02/2023
# - DESCRIPTION: Code to estimate the demographic resilience of Paramuricea 
#                clavata populations
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 
CodePath <- paste0(path,"/Code")

# Libraries

library(dplyr)
library(tidyverse)
library(popdemo)
library(popbio)
library(Matrix)

# Load the data

load(paste0(DataPath, "/Matrices.RData"))

# Load functions 

source(paste0(CodePath, "/DemoFunctions.R"))

# Join the different matrices 

mats_all <- mats_medes %>% dplyr::select(-vr) %>% 
  rbind(mats_pc %>% dplyr::select(-vr), 
        mats_m_old)

# Estimate the different demographic resilience metrics

dem_res <- mats_all %>% 
  # Transform into matrices
  mutate(matrix=map(matrix,
                    ~matrix(unlist(.x),
                            ncol=7, byrow=TRUE)),
         # Estimate its ergodicity, primitivity, and irreducibility
        # These properties are needed to trust the demographic resilience estimates
         ergodic=map(matrix,
                     ~isErgodic(.x)),
         primitive=map(matrix,
                       ~isPrimitive(.x)),
         irreducible=map(matrix,
                     ~isIrreducible(.x))) %>% 
  # Filter out those matrices that are ergodic, primitive, and irreducible
  filter(ergodic==TRUE & primitive==TRUE & irreducible==TRUE) %>% 
  # Remove the three years matrix in Port Cros
  filter(!(site=="Port-Cros"&year=="2001-2003")) %>% 
  # Calculate the demographic resilience metrics
  mutate(compensation= map(matrix,
                     ~maxamp(.x)),
         resistance=map(matrix,
                       ~reac(.x, bound = "lower")),
         speed_recovery=map(matrix,
                         ~dr(.x)),
         e_res=map(matrix, 
                   ~elas_res(.x, )),
         e_rec=map(matrix, 
                           ~elas_rec(.x)),
         # REMOVE BEFORE PUBLICATION
         time_recovery=map(matrix, 
                           ~return.time(.x)),
          r_up=map(matrix, ~reac(.x, bound = "upper")),
          r_low=map(matrix, ~reac(.x, bound = "lower"))) %>% 
  # unnest these estimates
  unnest(c(compensation, resistance, speed_recovery, time_recovery, 
           e_res, e_rec),names_sep = "_") 
  
# Save the data 

setwd(DataPath)
save(dem_res,file = "DemRes.RData")

