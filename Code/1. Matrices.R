# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Matrices.R         
# - DATE:        20/02/2023
# - DESCRIPTION: Code to build matrix population models from field data for
#                Paramuricea clavata populations
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(dplyr)
library(tidyverse)

# Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 
CodePath <- paste0(path,"/Code")

# Load data 

data_medes <- read.csv2(paste0(DataPath,"/Medes.csv"))
data_pc <- read.csv2(paste0(DataPath,"/PortCros.csv"))

# Load the function

source(paste0(CodePath, "/DemoFunctions.R"))

# Calculate vital rates ########################################################

# Predetermine the breaks, size classes, and fertility

breaks <- c(0, 0.3, 3, 10, 20, 30, 40, 400)
class.names <- c("0-0.3", "0.3-3", "3-10", "11-20", "21-30", "31-40", ">40")
fert <- c(0, 0, 0.29*0.7, 0.55*0.537, 0.95*0.781, 1*0.864, 1*0.900)

## Medes -----------------------------------------------------------------------

# We consider that those that have a 100% necrosis are dead

data_medes <- data_medes %>% 
  mutate(# Create a column with total mortalities
    T16=E16+N16, 
    T17=E17 + N17,
    T18=E18 + N18,
    T19=E19 + N19,
    T20=E20 + N20,
    T21=E21 + N21,
    T22=E22 + N22,
    # Assign 0 height to those with 100 mortality
    HC16=ifelse(is.na(T16), 0, ifelse(T16 == 100, 0, H16)),
    HC17=ifelse(is.na(T17), 0, ifelse(T17 == 100, 0, H17)),
    HC18=ifelse(is.na(T18), 0, ifelse(T18 == 100, 0, H18)),
    HC19=ifelse(is.na(T19), 0, ifelse(T19 == 100, 0, H19)),
    HC20=ifelse(is.na(T20), 0, ifelse(T20 == 100, 0, H20)),
    HC21=ifelse(is.na(T21), 0, ifelse(T21 == 100, 0, H21)),
    HC22=ifelse(is.na(T22), 0, ifelse(T22 == 100, 0, H22))) %>% 
  # Select only data from 2016 to 2021
  select(Estacio, HC16:HC22) %>% 
  filter(Estacio != "Medallot")

# Count number of colonies

data_medes %>% 
  pivot_longer(cols = 2:8,names_to = "year", values_to = "num") %>% 
  group_by(year) %>% 
  summarise(n=n())
  
# Estimate vital rates 
# First we exclude a few stations

v_rates1 <- data_medes %>%
  filter(Estacio != "Tasco Petit",
         Estacio != "Carall Bernat",
         Estacio != "La Vaca") %>%
  # Add a variable indicating the site
  mutate(site="Medes") %>% 
  # Nest our height values inside site to be able to apply functions
  nest(df =-site) %>% 
  # Use the function vital rates for each estacio and year transition
  mutate(y1617=map(df,
                   ~vital_rates(size=.$HC16, 
                           size_next = .$HC17,
                           size_next2 = .$HC18,
                           breaks = breaks,
                           class=class.names,
                           fert=fert, as_list=F))) %>% 
  # Remove the column df
  select(-df) %>%
  # Pivot longer the vital rates estimated per year (more rows)
  pivot_longer(cols = 2, names_to = "year", 
               values_to = "vital_rates")

# Now the rest 

v_rates2 <- data_medes %>%
  # Add a variable indicating the site
  mutate(site="Medes") %>% 
  # Nest our height values inside site to be able to apply functions
  nest(df =-site) %>% 
  # Use the function vital rates for each estacio and year transition
  mutate(y1718=map(df,
                   ~vital_rates(size=.$HC17, 
                                size_next = .$HC18,
                                size_next2 = .$HC19,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y1819=map(df,
                   ~vital_rates(size=.$HC18, 
                                size_next = .$HC19,
                                size_next2 = .$HC20,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y1920=map(df,
                   ~vital_rates(size=.$HC19, 
                                size_next = .$HC20,
                                size_next2 = .$HC21,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y2021=map(df,
                   ~vital_rates(size=.$HC20, 
                                size_next = .$HC21,
                                size_next2 = .$HC22,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y2122=map(df,
                   ~vital_rates(size=.$HC21, 
                                size_next = .$HC22,
                                size_next2 = NA,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F))) %>% 
# Remove the column df
select(-df) %>%
  # Pivot longer the vital rates estimated per year (more rows)
  pivot_longer(cols = 2:6, names_to = "year", 
               values_to = "vital_rates")

# Because during the first year there was only one recruit, the values are biased
# so we will add data from previous years. 
# First, we calculate the mean value of the vital rates
m_value <- v_rates2 %>% 
  # Unnest the vital rates to add two extra columns, vital_rate and value  
  unnest(vital_rates) %>% 
  # Group the variables according to vital rate
  group_by(vital_rate) %>% 
  # Estimate the mean value across the data set
  summarise(m=mean(value)) %>% 
  # Pivot the dataset to make it more accessible
  pivot_wider(names_from = vital_rate, values_from = m)

# We modify the vital rates of the first year for the recruit's survival and growth

v_rates1 <- v_rates1 %>% 
  unnest(vital_rates) %>% 
  mutate(value=ifelse(vital_rate=="g2", m_value$g2,
                      ifelse(vital_rate=="s1", m_value$s2,
                             ifelse(vital_rate=="s2", m_value$s2, value)))) %>% 
  nest(vital_rates=c(vital_rate, value))

# Join them

v_rates <- v_rates1 %>%
    rbind(v_rates2)

## Port-Cros -------------------------------------------------------------------

# Count number of colonies

data_pc %>% 
  select(height, height.1, height.2, height.3, height.4) %>% 
  rowwise() %>% 
  mutate(total=sum(height, height.1, height.2, height.3,height.4, na.rm = T)) %>% 
  filter(total!=0)

data_pc %>% 
  select(height.5, height.6, height.7, height.8, height.9) %>% 
  rowwise() %>% 
  mutate(total=sum(height.5, height.6, height.7, height.8, height.9, na.rm = T)) %>% 
  filter(total!=0)


# Estimate vital rates from each three 

v_rates_2 <- data_pc %>%
  # Rename all columns to lower case
  rename_all(tolower) %>% 
  # Remove the places where we don't have value for "lloc"
  drop_na(lloc) %>% 
  # Add a variable indicating the site
  mutate(site="Port-Cros") %>% 
  # Nest our height values inside site to be able to apply functions
  nest(df =-site) %>% 
  # Use the function vital rates for each lloc and year transition
  mutate(y9900=map(df,
                   ~vital_rates(size=.$height, 
                                size_next = .$height.1,
                                size_next2 = .$height.2,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0001=map(df,
                   ~vital_rates(size=.$height.1, 
                                size_next = .$height.2,
                                size_next2 = .$height.3,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0102=map(df,
                   ~vital_rates(size=.$height.2, 
                                size_next = .$height.3,
                                size_next2 = NA,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0506=map(df,
                   ~vital_rates(size=.$height.5, 
                                size_next = .$height.6,
                                size_next2 = .$height.7,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0607=map(df,
                   ~vital_rates(size=.$height.6, 
                                size_next = .$height.7,
                                size_next2 = .$height.8,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0708=map(df,
                   ~vital_rates(size=.$height.7, 
                                size_next = .$height.8,
                                size_next2 = .$height.9,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)),
         y0809=map(df,
                   ~vital_rates(size=.$height.8, 
                                size_next = .$height.9,
                                size_next2 = NA,
                                breaks = breaks,
                                class=class.names,
                                fert=fert, as_list=F)))%>% 
  # Remove the variable df 
  select(-df) %>% 
  # Pivot longer the vital rates estimated per year (more rows)
  pivot_longer(cols = 2:8, 
               names_to = "year", 
               values_to = "vital_rates")


# Matrix building ##############################################################

## Create a baseline matrix ----------------------------------------------------

Asim <- expression(
  0,           0,           0,                f4,                 f5,                   f6,                   f7,
  s1*g1,       s2*(1-g2),   0,                 0,                  0,                    0,                    0,
  0,           s2*g2,       s3*(1-g3),         s4*(1-g4)*h4,       s5*(1-g5)*h5*hl5,     0,                    0,
  0,           0,           s3*g3,             s4*(1-g4)*(1-h4),   s5*(1-g5)*h5*(1-hl5), s6*(1-g6)*h6*hl6,     0,
  0,           0,           0,                 s4*g4,              s5*(1-g5)*(1-h5),     s6*(1-g6)*h6*(1-hl6), s7*h7*hl7,
  0,           0,           0,                 0,                  s5*g5,                s6*(1-g6)*(1-h6),     s7*h7*(1-hl7),
  0,           0,           0,                 0,                  0,                    s6*g6,                s7*(1-h7))

## Medes -----------------------------------------------------------------------

mats_medes <- v_rates %>% 
  # Unnest the vital rates to have them as two separate colums
  unnest(vital_rates) %>% 
  # Pivot the data frame to put the vital rates as columns
  pivot_wider(names_from = vital_rate,
              values_from = value) %>% 
  # Replace all the infinite values with NA
  mutate_if(is.numeric, ~ifelse(is.infinite(.), NA, .)) %>%  
  # Replace all the NA and missing values with the mean value of all other sites
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) %>%
  # Pivot back to longer format (more rows less columns)
  pivot_longer(cols = 3:34, 
               names_to = "vital_rate", 
               values_to = "value")

# Now we calculate the mean value per vital rate except 2021 

mean_v <- mats_medes %>% 
  filter(year != "y2021") %>% 
  group_by(vital_rate) %>% 
  summarize(mean_v = mean(value, na.rm=T)) %>% 
  filter(str_detect(vital_rate, "f")) %>% 
  mutate(year="y2021")

# Now we finish the data frame 

mats_medes <- mats_medes %>% 
  # Substitute values of f3, f4, f5, f6 and f7
  # We join the mean values
  left_join(mean_v) %>% 
  # We add the value from the mean_v when available
  mutate(value = ifelse(!is.na(mean_v), mean_v, value)) %>% 
  # Remove mean_v
  select(-mean_v) %>% 
  # Nest the vital rates inside year
  nest(vr =-c(site, year)) %>% 
  # Split the columns of the names of the vital rates and their value into a list
  mutate(vr=map(vr, function(vr) split(vr, f = vr$vital_rate) %>% 
                  map("value")),
         # Use the list of vital rates to create a matrix from the symbolic matrix
         matrix=map(vr,
                    ~matrix(sapply(Asim, eval, .),
                            ncol=7, byrow=TRUE,
                            dimnames=list(class.names, 
                                          class.names)))) %>% 
  # Modify year 
  mutate(year=ifelse(year=="y1617", "2016-2017",
                     ifelse(year=="y1718", "2017-2018",
                            ifelse(year=="y1819", "2018-2019",
                                   ifelse(year=="y1920", "2019-2020", 
                                          ifelse(year=="y2021", "2020-2021", 
                                                 ifelse(year=="y2122", "2021-2022", 
                                                        year)))))))



## Port-Cros -------------------------------------------------------------------

mats_pc <- v_rates_2 %>% 
  # Unnest the vital rates to have them as two separate colums
  unnest(vital_rates) %>% 
  # Pivot the data frame to put the vital rates as columns
  pivot_wider(names_from = vital_rate,
              values_from = value) %>% 
  # Replace all the infinite values with NA
  mutate_if(is.numeric, ~ifelse(is.infinite(.), NA, .)) %>% 
    # Replace all the NA and missing values with the mean value of all other sites
  mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE))) %>%
  pivot_longer(cols = 3:37, 
               names_to = "vital_rate", 
               values_to = "value") %>% 
  # Separate the duplicated values
  unnest(value) %>%  
  # Group by site year and vital_rate
  group_by(site, year, vital_rate) %>%
  # Keep the second repeated value
  slice_tail(n = 1)%>% 
  # Ungroup
  ungroup()
  
# Now we calculate the mean value per vital rate except 2021 

mean_vals <- mats_pc %>%
  filter(str_detect(vital_rate, "f")) %>%
  group_by(vital_rate) %>%
  summarize(mean_value = mean(value))

# Now we finish the data frame 

mats_pc <- mats_pc %>% 
  # Substitute values of f3, f4, f5, f6 and f7 when they are 0
  # We join the mean values
  left_join(mean_vals, by="vital_rate") %>% 
  # We add the value from the mean_values when the present value is 0
  mutate(value = ifelse(vital_rate %in% c("f3", "f4", "f5", "f6", "f7") & value == 0,
                        mean_value, value)) %>%
  # Remove mean_v
  select(-mean_value) %>%
  # Nest the vital rates inside year and lloc
  nest(vr =-c(site, year)) %>% 
  # Split the columns of the names of the vital rates and their value into a list
  mutate(vr=map(vr, function(vr) split(vr, f = vr$vital_rate) %>% 
                  map("value")),
         # Use the list of vital rates to create a matrix from the symbolic matrix
         matrix=map(vr,
                    ~matrix(sapply(Asim, eval, .),
                            ncol=7, byrow=TRUE,
                            dimnames=list(class.names, 
                                          class.names)))) %>% 
  # Modify year 
  mutate(year=ifelse(year=="y9900", "1999-2000",
                     ifelse(year=="y0001", "2000-2001",
                            ifelse(year=="y0102", "2001-2003", 
                                   ifelse(year== "y0506", "2005-2006",
                                          ifelse(year=="y0607", "2006-2007", 
                                                 ifelse(year=="y0708", "2007-2008",
                                                        ifelse(year=="y0809", 
                                                               "2008-2009", year))))))))

## Medes 2001-2004 from Linares et al 2007 -------------------------------------

# 2001-2002

a0102 <- matrix(c(0,     0,     0, 0.003, 0.089, 0.378, 0.810,
                  0.667, 0.380, 0.033, 0,     0,     0,     0,
                  0, 0.292, 0.730, 0.043,     0,     0,     0,
                  0,     0, 0.120, 0.780, 0.162, 0.060,     0,
                  0,     0,     0, 0.080, 0.736, 0.092, 0.038,
                  0,     0,     0,     0, 0.060, 0.660, 0.057,
                  0,     0,     0,     0,     0, 0.092, 0.740),
                ncol=7, byrow=TRUE,
                dimnames=list(class.names, class.names))

# 2002-2003

a0203 <- matrix(c(0,     0,     0, 0.003, 0.089, 0.378, 0.810,
                  0.636, 0.500,  0,    0,     0,     0,     0,
                  0, 0.136, 0.860, 0.020,     0,     0,     0,
                  0,     0, 0.102, 0.850, 0.003,     0,     0,
                  0,     0,     0, 0.101, 0.750, 0.120, 0.000,
                  0,     0,     0,     0, 0.078, 0.850, 0.100,
                  0,     0,     0,     0,     0, 0.066, 0.880),
                ncol=7, byrow=TRUE,
                dimnames=list(class.names, class.names))     

# 2003-2004 

a0304 <- matrix(c(0,     0,     0, 0.003, 0.089, 0.378, 0.810,
                  0.750, 0.570, 0.020, 0,     0,     0,     0,
                  0, 0.179, 0.850, 0.040,     0,     0,     0,
                  0,     0, 0.040, 0.870, 0.007,     0,     0,
                  0,     0,     0, 0.034, 0.790, 0.160, 0.020,
                  0,     0,     0,     0, 0.017, 0.683, 0.020,
                  0,     0,     0,     0,     0, 0.032, 0.930),
                ncol=7, byrow=TRUE,
                dimnames=list(class.names, class.names))

# Create a data frame with the years and the location

mats_m_old <- data.frame(site="Medes old", 
                         year=c("2001-2002", 
                                "2002-2003", 
                                "2003-2004"))

# Combine the matrices into a single data frame using rbind() and cbind()

mats_m_old <- cbind(mats_m_old, 
                    rbind(data.frame(matrices=I(list(a0102))), 
                          data.frame(matrices=I(list(a0203))), 
                          data.frame(matrices=I(list(a0304)))))

# Nest the data set 

mats_m_old <- mats_m_old %>% 
  nest(matrix = matrices)

# Classify now Port Cros as old or new 

mats_pc <- mats_pc %>% 
  mutate(site=ifelse(year%in%c("1999-2000", 
                               "2000-2001"),
                     "Port-Cros old", "Port-Cros"))


# Save the data ----------------------------------------------------------------

setwd(DataPath)
save(list = ls(pattern = 'mats_'), 
     file = "Matrices.RData")
  