# --------------------------------------------------------------------------------------- #
# - FILE NAME:   LTRE.R         
# - DATE:        12/04/2023
# - DESCRIPTION: Code to calculate LTRE of Paramuricea clavata populations 
#                exposed and not exposed to MHWs. 
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
library(lefko3)

#Personal theme for ggplot

theme_set(theme(panel.background=element_blank(), 
                strip.background=element_blank(), 
                axis.line=element_line("black"), 
                axis.ticks=element_line("black"), 
                axis.text=element_text(colour="black", size=12), 
                axis.title=element_text(size=16),
                panel.grid=element_blank(),
                strip.text=element_text(size=14),
                legend.key=element_blank(), 
                axis.line.x=element_line("black"),
                axis.line.y=element_line("black"),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                plot.margin = margin(0.1, 0, 0, 0, "cm"),
                plot.tag = element_text(face = 'bold', size = 16)))    

# Load the data

load(paste0(DataPath, "/Matrices.RData"))
load(paste0(DataPath, "/temperature.RData"))

# Join all the data sets to make the calculations easier 

mats_all <- rbind(mats_m_old %>% 
                    unnest(matrix) %>% 
                    rename(matrix=matrices) %>% 
                    mutate(vr=NA),
                  mats_medes,
                  mats_pc) %>% 
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
  filter(!(site=="Port-Cros"&year=="2001-2003"))

# Join the different MHWs categories

categories <- rbind(mhw_cat %>% mutate(site="Port-Cros"), 
                    mhw_cat_m%>% mutate(site="Medes")) 

# Join with the demographic data

mats <- mats_all %>% 
  rename(location=site,
         Year=year) %>% 
  mutate(site=gsub(" .*$", "", location),
         year=gsub("-.*", "",Year))  %>% 
  left_join(categories %>% 
              mutate(year=format(peak_date, "%Y")) %>% 
              dplyr::select(site, year, category, i_max)) %>% 
  mutate(category=ifelse(is.na(category), 
                         "None", category),
         i_max=ifelse(is.na(i_max), 0, i_max)) %>% 
  group_by(site, year) %>% 
  filter(i_max==max(i_max, na.rm = T)) %>% 
  mutate(cat= ifelse(category%in%c("None", "I Moderate"), "No-MHW", "MHW"))

# We will define three different sets of matrices that will be used for the following
# analyses

ref_mats <- mats %>% filter(cat=="No-MHW")
com_mats <- mats %>% filter(cat=="MHW")

# Small Noise Approximation for LTRE MHW/No-MHW --------------------------------

# Create the size vector with the sizes of the stages, for this we use the midpoints

sizevector <- c(0.15, (0.3+3)/2,(3+10)/2,(11+20)/2,(21+30)/2,(31+40)/2, 40)

# Create the stage vector with the names

stagevector <- rownames(mats_pc$matrix[[1]])

# Reproductive status, binary value on whether the stage is reproductive or not 

repvector <- as.integer(mats_pc$matrix[[1]][1,] != 0)

# Observed status, binary value on whether the stage is observed or not

obsvector <- c(0,1,1,1,1,1,1)

# Maturity status, binary value on whether the stage is morphologically mature or not

matvector <- c(0,1,1,1,1,1,1)

# Propagule status, binary value on whether the stage is a propagule

propvector <- c(1,0,0,0,0,0,0)

# Size of the bins

binsize <- c(0.15,1.35,3.5,4.5,4.5,4.5,0)

# Create the population data frame

pop_frame <- sf_create(sizes = sizevector,stagenames = stagevector,
                       repstatus = repvector,obsstatus = obsvector,
                       matstatus = matvector, propstatus = propvector,
                       binhalfwidth = binsize)

# Reference matrices

ref_mats <- create_lM(mats = ref_mats$matrix,
                      stageframe = pop_frame,
                      historical = FALSE,
                      yearorder = ref_mats$year)

# MHW matrices

com_mats <- create_lM(mats = com_mats$matrix,
                      stageframe = pop_frame,
                      historical = FALSE,
                      yearorder = com_mats$year)


# Small Noise Approximation LTRE

ltre_res <- ltre3(mats = com_mats,
                  refmats = ref_mats,
                  sna_ltre = T)

# Lambda difference

exp(ltre_res$r_values_m)
exp(ltre_res$r_value_ref)

# Summarise the results

ltre_summary <- summary(ltre_res)

# Create the data frame with the contributions of the different demographic processes
# to the mean matrix

ltre_dat_mean <- data.frame(cont_mean=rbind(cbind(ltre_summary$ahist_mean$matrix1_pos),
                                            cbind(ltre_summary$ahist_mean$matrix1_neg)),
                            demogr=ltre_summary$ahist_mean$category)

# Create a data frame with the absolute contributions of the mean, elasticity, 
# covariartion and correlation

ltre_dat_abs <- data.frame(cont=rbind(ltre_summary$overall$means_positive,
                                      ltre_summary$overall$elas_positive,
                                      ltre_summary$overall$cv_positive,
                                      ltre_summary$overall$cor_positive,
                                      ltre_summary$overall$means_negative,
                                      ltre_summary$overall$elas_negative,
                                      ltre_summary$overall$cv_negative,
                                      ltre_summary$overall$cor_negative),
                           component=c("means", "elasticity", 
                                       "CV", "correlation"))

# Save the results
setwd(ResultPath)
save(ltre_dat_abs, ltre_dat_mean,file = "LTRE.RData")
