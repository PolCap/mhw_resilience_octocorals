# --------------------------------------------------------------------------------------- #
# - FILE NAME:   ResilienceTemperature.R         
# - DATE:        25/05/2023
# - DESCRIPTION: Code to study the influence of temperature on the resilience of 
#                Paramuricea clavata populations
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
library(Rage)

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
                # text = element_text(family = "futura"),
                plot.tag = element_text(face = 'bold', size = 16)))  

# Load the data

load(paste0(DataPath, "/DemRes.RData"))
load(paste0(DataPath, "/temperature.RData"))

# Temporal trends --------------------------------------------------------------

# Gradual temperatures

# Select the climatology from Port-Cros, add the name of the site

temp <- mhw_pc$climatology %>% 
  dplyr::select(t,temp) %>% 
  mutate(site="Port-Cros")

# Estimate the mean

temp <- temp %>% 
  mutate(year=format(t, "%Y")) %>% 
  group_by(site, year) %>% 
  summarise(temp=mean(temp, na.rm=T)) %>% 
  rbind(temp_m2 %>% 
          mutate(year=as.character(year),
                 site="Medes") %>% 
          dplyr::select(site, year, temp))

# Join with resilience metrics

dem_grad <- dem_res %>% 
  rename(location=site,Year=year) %>% 
  mutate(site=gsub(" .*$", "", location),
         year=gsub(".*\\-", "",Year)) %>% 
  left_join(temp, by=c("site", "year"))

# We calculate the duration and intensity mhws by years and site, then we 
# join it

events <- rbind(mhw_m$event %>% 
                  mutate(year = format(date_start, "%Y")) %>% 
                  group_by(year) %>% 
                  summarise(duration = sum(duration),
                            intensity = max(intensity_max_abs)) %>% 
                  mutate(site= "Medes"), 
                mhw_pc$event %>% 
                  mutate(year = format(date_start, "%Y")) %>% 
                  group_by(year) %>% 
                  summarise(duration = sum(duration),
                            intensity = max(intensity_max_abs)) %>% 
                  mutate(site= "Port-Cros"))

# Modify site and year variables in the resilience dataframe 
# and join with events

dem <- dem_res %>% 
  rename(location=site,Year=year) %>% 
  mutate(site=gsub(" .*$", "", location),
         year=gsub("-.*", "",Year)) %>% 
  left_join(events) %>% 
  mutate(intensity=ifelse(is.na(intensity), 0, intensity),
         duration=ifelse(is.na(duration), 0, duration))

# Heatwaves categories ---------------------------------------------------------

# Join the different MHWs categories

categories <- rbind(mhw_cat %>% mutate(site="Port-Cros"), 
                    mhw_cat_m%>% mutate(site="Medes")) 

# Join with the demographic data

dem_cat <- dem_res %>% 
  rename(location=site,Year=year) %>% 
  mutate(site=gsub(" .*$", "", location),
         year=gsub(".*\\-", "",Year))  %>% 
  left_join(categories %>% 
              mutate(year=format(peak_date, "%Y")) %>% 
              dplyr::select(site, year, category, i_max, duration)) %>% 
  mutate(category=ifelse(is.na(category), 
                         "None", category)) %>% 
  # distinct(site, year, .keep_all = T) %>% 
  dplyr::select(site, year, category, resistance, compensation, 
         time_recovery, speed_recovery, 
         e_res_fecundity, e_res_growth,e_res_shrinkage,e_res_stasis, 
         e_rec_fecundity, e_rec_growth,e_rec_shrinkage,e_rec_stasis,
         i_max, duration) %>% 
  mutate(intensity=ifelse(is.na(i_max), 0, i_max),
         duration=ifelse(is.na(duration), 0, duration)) %>% 
  group_by(site, year) %>% 
  filter(duration==max(duration, na.rm = T)) %>% 
  mutate(cat= ifelse(category%in%c("None", "I Moderate"), "No-MHW", "MHW"),
         cat=factor(cat, levels=c("No-MHW", "MHW")),
         category=factor(category,
                                c("None", "I Moderate","II Strong", "IV Extreme")))  
  

# Save data sets for posterior modelling

setwd(ResultPath)
save(dem_cat, dem_grad, file = "ResTemp.RData")

