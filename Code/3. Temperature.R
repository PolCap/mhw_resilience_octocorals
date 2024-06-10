# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Temperature.R         
# - DATE:        12/04/2023
# - DESCRIPTION: Code to calculate MHWs from temperature data. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Working directories

setwd(gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path)))

# Libraries

library(dplyr)
library(tidyverse)
library(patchwork)
library(heatwaveR)

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
                plot.margin = margin(0.1, 0, 0, 0, "cm")))  

# Load the data

temp_pc <- read.csv2("Data/T_Port_Cros.csv")
temp_m <- read.delim("Data/T_Medes.txt",check.names = F)

# Port-Cros --------------------------------------------------------------------

# Transform the data

temp_pc <- temp_pc %>% 
  mutate(t=as.Date(Date,"%d/%m/%y")) %>% 
  rename(temp=X25.m) %>% 
  select(t, temp)

# Create the climatology time series

ts <- ts2clm(temp_pc, climatologyPeriod = c("1999-06-17", "2009-07-24"))

# Filter the summer days

ts <- ts %>% 
  filter(doy >= 153 & doy <= 335)

# Detect the events in a time series

mhw_pc <- detect_event(ts)

# Detect categories 

mhw_cat <- category(mhw_pc, S = FALSE, name = "PC")

# Medes ------------------------------------------------------------------------

# Transform the data

temp_m <- temp_m %>% 
  pivot_longer(cols = 3:10, names_to = "depth", 
               values_to = "temp") %>% 
  filter(depth==20) %>% 
  mutate(t=as.Date(Date,"%d/%m/%Y")) %>%
  group_by(t) %>% 
  summarise(temp=mean(temp, na.rm=T)) %>% 
  select(t, temp)

# Create the climatology time series

ts_m <- ts2clm(temp_m, climatologyPeriod = c("2002-07-18", "2012-10-19"))

# Filter the summer days

ts_m <- ts_m %>% 
  filter(doy >= 153 & doy <= 335)

# Detect the events in a time series

mhw_m <- detect_event(ts_m)

# Detect categories 

mhw_cat_m <- category(mhw_m, S = FALSE, name = "Medes")

# Save the data 

setwd(DataPath)
save(file = "temperature.RData", mhw_m, mhw_pc, mhw_cat, mhw_cat_m)
