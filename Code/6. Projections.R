# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Projections.R         
# - DATE:        28/06/2023
# - DESCRIPTION: Code to project the P. clavata populations in different MHWs 
#                scenarios.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Working directories

setwd(gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path)))

# Libraries

library(dplyr)
library(tidyverse)
library(popdemo)
library(popbio)
library(Matrix)
library(patchwork)
library(ggdist)
library(ggrepel)
library(foreach)
library(doParallel)

# Personal theme for ggplot

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

load("Data/ResTemp.RData")

# Join the two datasets

dem <- dem_cat %>% left_join(dem_grad)

# Load functions 

source("Code/DemoFunctions.R")

# Parameters for the simulations

probs <- seq(0, 1, by=0.2)# Probabilities of disturbance
tmax <- 100 # Time of projection
reps <- 1000 # Number of simulations
n0 <- stable.stage(dem$matrix[11][[1]])*100 # Initial population

# Set the number of cores to utilize in parallel

num_cores <- 12

# Create a parallel backend using doParallel

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Perform the loop in parallel

results <- foreach(i = probs, .combine = rbind,
                   .packages = c("tidyverse", "popdemo")) %dopar% {
                     n <- pop_proj(mats = dem$matrix[dem$cat == "No-MHW"],
                                   mats_mhw = dem$matrix[dem$cat == "MHW"],
                                   prob = i,
                                   n0 = n0,
                                   tmax = tmax,
                                   reps = reps)
                     }

# Stop the parallel backend

stopCluster(cl)

# Save the results

setwd(ResultPath)
save(results, file = "Simulations.RData")

# Figure Population abundance ----------------------------------------------

# Find the first time that the population abundance goes below 1

(position <- results %>% 
    mutate(probs=paste("Every", prob*100, "years"),
           probs=fct_reorder(probs, prob)) %>%
    group_by(probs, rep, time) %>% 
    mutate(above.qe= n >= 10,
           alife=cumprod(above.qe)) %>% 
    group_by(probs, time) %>% 
    summarise(prob= 1-(sum(alife, na.rm = T)/reps)) %>% 
    filter(prob==1) %>% 
    group_by(probs) %>% 
    summarise(time=min(time)))

# Create a gradient of colours 

pal <-  colorRampPalette(c("#D9AA1E", "#A6122D"))
cols <- pal(length(probs))

# Plot it

(p3 <- results %>%
   filter(rep<250) %>% 
    mutate(probs=paste("Every", prob*100, "years"),
           probs=fct_reorder(probs, prob)) %>%
    ggplot(aes(x=time, y=n,
              colour=probs))+
    geom_line(aes(group=interaction(probs, rep)),
              alpha=0.1) +
    geom_hline(aes(colour=probs, 
               yintercept = 1), 
               linetype="dashed")+
    facet_wrap(~probs)+
    geom_text_repel(data = position,
                    aes(x=time, y= 1,
                        label=paste(time, "years")),
                    segment.color = 'black',
                    box.padding = 2.5,
                    nudge_x = .15,
                    nudge_y = 1, 
                    size=5) +
    scale_colour_manual(values=cols)+
    labs(y= "Abundance",
         x="Time")+
    theme(legend.position = "none"))

