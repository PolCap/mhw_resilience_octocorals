# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupplementaryFigures.R         
# - DATE:        01/05/2024
# - DESCRIPTION: Code to create the suppplementary figures.
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila@ub.edu)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(dplyr)
library(tidyverse)
library(sf)
library(marmap)
library(sp)
library(maptools)
library(rgeos)
library(maps)
library(rgdal)
library(raster)
library(gridBase)
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

# Working directory

setwd(gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path)))

# Figure S1: Map ---------------------------------------------------------------

# Sampling sites 

places <- data.frame(site_name=c("Pedra de Déu", "Pota del Llop", 
                                 "la Vaca Nord", "Tascó Gros", "Tascó Petit", 
                                 "Carall Bernat", "Gabinière", "Montremian"), 
                     latitude=c(42.050137, 42.049099, 
                                42.047581, 42.041961, 42.041352, 
                                42.041518, 42.987237, 43.018625), 
                     longitude=c(3.224590, 3.226176, 
                                 3.225539, 3.226916, 3.226649,
                                 3.228496, 6.393867, 6.362907))


# Transform data points into geographic objects

Sites_geo <- places %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Get bathymetry data

bathy <- getNOAA.bathy(-10,20, 33, 46.5, res = 1, keep = TRUE)
ggbathy <- fortify(bathy)

# Get countries outline

pays <- rnaturalearth::ne_countries(
  country = c("France", "Spain", "Italy","Andorra", "Portugal",
              "Libya", "Bosnia and Herzegovina","Slovenia","Greece",
              "Albania","Montenegro","Republic of Serbia", "Kosovo",
              "Hungary",
              "Morocco","Algeria", "Tunisia", "Croatia", "Switzerland"),
  scale = "large", returnclass = "sf")

#Our transformation function

scaleFUN <- function(x) sprintf("%.2f", x)

# Main map ---------------------------------------------------------------------

# Write the locations 

names <- data.frame(Lon=c(2.9, 5.9), 
                    Lat=c(42.2, 43.2), 
                    label=c("b", "c"))
# Make the plot

(pl <- ggplot(data = pays) +
    geom_contour(data = ggbathy, aes(x = x, y = y, z = z),
                 binwidth = 800, color = "grey80", size = 0.1) +
    geom_sf() +
    geom_rect(
      xmin = 2.9,
      ymin = 41.8,
      xmax = 3.5,
      ymax = 42.2,
      fill = NA, 
      colour = "black",
      size = 0.6)+
    geom_rect(
      xmin = 5.9,
      ymin = 42.8,
      xmax = 6.5,
      ymax = 43.2, 
      fill = NA, 
      colour = "black",
      size = 0.6)+
    coord_sf(xlim = c(-8,19), 
             ylim = c(34, 46)) +
    geom_text(data=names, aes(x=Lon, y=Lat,
                              label=label), 
              size = 8, hjust=1.1, vjust=1.1) +  
    ggspatial::annotation_scale(
      location = "tr",
      bar_cols = c("grey60", "white")) +
    ggspatial::annotation_north_arrow(
      location = "tl", which_north = "true",
      pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm"))+
    labs(x="", y="")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# Inset 1 ----------------------------------------------------------------------

# Load data 
# Coastline 

coast <-  readShapePoly(paste0(CodePath, "/maps/shp/Costa.shp"),
                        proj4string = CRS("+proj=utm +zone=31 +datum=WGS84"))
ggcoast <- st_as_sf(coast)

ggcoast <- st_crop(ggcoast,
                   xmin = 518005,
                   ymin = 4654195,
                   xmax = 519050,
                   ymax = 4655600)

# Bathymetry

bati <- readShapeLines(paste0(CodePath, "/maps/shp/Bati.shp"),
                       proj4string = CRS("+proj=utm +zone=31 +datum=WGS84"))

ggbat <- st_as_sf(bati)

# Filter the points

p <- places %>% filter(site_name%in%c("Medallot", "Pedra de Déu", "Pota del Llop", 
                                      "la Vaca Nord", "Tascó Gros", "Tascó Petit", 
                                      "Carall Bernat"))

# Create a SpatialPoints object with your decimal coordinates

coords <- SpatialPoints(matrix(c(p$longitude, p$latitude), ncol = 2), 
                        proj4string = CRS("+proj=longlat +datum=WGS84"))

# Transform the coordinates to the desired UTM projection.

utm_proj <- CRS("+proj=utm +zone=31 +datum=WGS84")  # Assuming UTM Zone 18 for New York City
coords_utm <- spTransform(coords, utm_proj)
points <- st_as_sf(coords_utm)

# Inset 

(p1 <- ggplot() +
    # geom_sf(data = ggbat, color = "#00000030", size=.1)+
    geom_sf(data = ggcoast, 
            color = "black") +
    ggspatial::annotation_scale(
      location = "tr",
      bar_cols = c("grey60", "white"))+ 
    geom_sf(data = points, shape = 21, fill="steelblue", 
            size=3, alpha=.9) +
    labs(x="", y="")+
    coord_sf(xlim =c(518005, 519050), 
             ylim = c(4654195, 4655600)) +
    scale_y_continuous(breaks = c(42.040, 42.045, 42.050)) + 
    scale_x_continuous(breaks = c(3.220,3.225, 3.230)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# Inset 2 ----------------------------------------------------------------------

# Load european coastline 

europe <- read_sf(paste0(CodePath,'/maps/shp/Europe_coastline_poly.shp'))

# Transform to sf

eu_coast <- st_transform(europe, crs = 4326)

# Crop the european coast to fit Port-Cros

portcros_coast <- st_crop(eu_coast,
                          xmin = 6.30,
                          ymin = 42.95,
                          xmax = 6.65,
                          ymax = 43.07)
# Subset the points

p <- places %>% filter(site_name%in%
                         c("Gabinière", "Montremian")) 

# Create a SpatialPoints object with your decimal coordinates

coords <- SpatialPoints(matrix(c(p$longitude,
                                 p$latitude), ncol = 2), 
                        proj4string = CRS("+proj=longlat +datum=WGS84"))

# Transform the coordinates to the desired UTM projection.

utm_proj <- CRS("+proj=utm +zone=31 +datum=WGS84")  # Assuming UTM Zone 18 for New York City
coords_utm <- spTransform(coords, utm_proj)
points <- st_as_sf(coords_utm)

# Plot 

(p2 <- ggplot() +
    geom_sf(data = portcros_coast, 
            color = "black") +
    ggspatial::annotation_scale(
      location = "br",
      bar_cols = c("grey60", "white"))+ 
    geom_sf(data = points,
            shape = 21,
            fill="steelblue",
            size=3, alpha=.9) +
    labs(x="", y="")+
    coord_sf(xlim =c(6.34, 6.52), 
             ylim = c(42.97, 43.06)) +
    scale_y_continuous(breaks = c(43.05, 43.00, 42.95)) + 
    scale_x_continuous(breaks = c(6.40, 6.50)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# Combine the maps -------------------------------------------------------------

(figS1 <- pl+(p1/p2)+
   plot_annotation(tag_levels = "a")&
   theme(plot.margin = unit(c(0,0,0,0), "cm")))

ggsave(figS1,
       filename = "Results/Figure S1.pdf",
       width = 14, height = 10)


# Figure S2: Elasticities ------------------------------------------------------
# Load the data

load("Results/Matrices.RData")

# Join the different matrices 

mats_all <- mats_medes %>% dplyr::select(-vr) %>% 
  rbind(mats_pc %>% dplyr::select(-vr), 
        mats_m_old)

# We will need to define the U and F matrices first, for that
# we create a function to create them from our matrices

mat_decomp <- function(x){
  matF <- matU <- x
  matF[-1,] <- 0
  matU[1,] <- 0
  return(tibble(matF, matU))
}

# Also we will create a function to summarise the different contributions decomposed
# by the different parts of the matrix, corresponding to different demographic
# processes.

mat_comp <- function(mat) {
  # decompose the matrix into matF and matU
  mats<- mat_decomp(mat)
  # Contributions to the reproduction
  Fecundity <- sum(mats$matF)
  # Contributions to stasis
  Stasis <- sum(diag(mats$matU))
  # Contributions to shrinkage
  Shrinkage <- sum(mats$matU[upper.tri(mats$matU)])
  # Contributions to growth
  Growth <- sum(mats$matU[lower.tri(mats$matU)])
  return(data.frame(Fecundity, Growth,
                    Shrinkage, Stasis))
}

# Calculate elasticities 

elas <- mats_all %>% 
  #Split the matrices using own function mat_decomp
  mutate(# Correct the matrices
    matrix=map(matrix,
                    ~matrix(unlist(.x),
                            ncol=7, byrow=TRUE)),
    # Matrix elasticities
    el_m=map(matrix,
                ~elasticity(.x)),
    
    # Decompose the elasticity matrices and add their elasticities
    el_elas =map(el_m,
                 ~mat_comp(.x)))

 # Plot it

(g1 <- elas %>%
  # Take the vital rate elasticity values into a wide format
  unnest_wider(el_elas) %>% 
  # Pivot longer to we have it in a ggplot friendly format
  pivot_longer(cols=c("Stasis", "Growth", "Shrinkage", "Fecundity"),
               names_to = "demo_proc",
               values_to = "elast") %>%
  # Modify some parameters for the plot 
  mutate(site=gsub(" old", "", site),
         year=gsub("\\-.*","", year)) %>% 
  # We plot it now
  ggplot(aes(x=demo_proc, y=elast, group=demo_proc, fill=demo_proc))+
    geom_boxplot()+ 
    labs(x="Demographic process", y="Elasticity value")+
    scale_fill_manual(name=NULL, values = c("#335c67","#fff3b0",
                                          "#e09f3e", "#9e2a2b"),
                      guide=NULL))

# Save the plot

ggsave("Results/Figure S2.pdf", g1,
       width = 6, height = 4)
  
# Figure S3: Characterise the thermal regime -----------------------------------
# Load the data

temp_pc <- read.csv2("Data/T_Port_Cros_25m.csv")
temp_m <- read.delim("Data/T_Medes.txt",check.names = F)

# Port-Cros --------------------------------------------------------------------

# Transform the data

temp_pc <- temp_pc %>% 
  mutate(t=as.Date(Date,"%d/%m/%y")) %>% 
  rename(temp=X25.m) %>% 
  dplyr::select(t, temp)

# Create the climatology time series

ts <- ts2clm(temp_pc, climatologyPeriod = c("1999-06-17", "2009-07-24"))

# Detect the events in a time series

mhw_pc <- detect_event(ts)

# Medes ------------------------------------------------------------------------

# Transform the data

temp_m <- temp_m %>% 
  pivot_longer(cols = 3:10, names_to = "depth", 
               values_to = "temp") %>% 
  filter(depth==20) %>% 
  mutate(t=as.Date(Date,"%d/%m/%Y")) %>%
  group_by(t) %>% 
  summarise(temp=mean(temp, na.rm=T)) %>% 
  dplyr::select(t, temp)

# Create the climatology time series

ts_m <- ts2clm(temp_m, climatologyPeriod = c("2002-07-18", "2012-10-19"))

# Detect the events in a time series

mhw_m <- detect_event(ts_m)

# Plot
# Medes 

(p1 <- ts_m %>% 
  left_join(mhw_m$event %>% 
              dplyr::select(date_peak, intensity_cumulative), 
            by = c("t"="date_peak")) %>% 
  ggplot(aes(x = t, y = intensity_cumulative, colour=intensity_cumulative)) +
  geom_lolli() +
  geom_vline(xintercept = as.Date("2002-01-01"))+
  geom_vline(xintercept = as.Date("2005-01-01"))+
  annotate(geom = "text", 
           colour = "black", 
           x = as.Date("2003-07-01"), 
           y = 170,
           label = "First period\nof the study",
           hjust = 0.5) +
  geom_vline(xintercept = as.Date("2016-01-01"))+
  geom_vline(xintercept = as.Date("2023-01-01"))+
  annotate(geom = "text", 
           colour = "black", 
           x = as.Date("2019-07-01"), 
           y = 170,
           label = "Second period \nof the study",
           hjust = 0.5) +
  scale_colour_gradient2("MHW cumulative \nintensity", 
                         low = "#ffc866",
                         mid = "#ff6900",
                         high = "#9e0000",
                         limits = c(0, 170))+
    labs(x = NULL, 
         y = expression(atop("Cumulative intensity", 
         paste("(", degree, "C" %*%"days)"))))+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,190))+
    theme(legend.title = element_text(hjust = .5))) 

# Port-Cros

(p2 <- ts %>% 
    left_join(mhw_pc$event %>% 
                dplyr::select(date_peak, intensity_cumulative), 
              by = c("t"="date_peak")) %>% 
    ggplot(aes(x = t, y = intensity_cumulative, colour=intensity_cumulative)) +
    geom_lolli() +
    geom_vline(xintercept = as.Date("1999-01-01"))+
    geom_vline(xintercept = as.Date("2004-01-01"))+
    annotate(geom = "text", 
             colour = "black", 
             x = as.Date("2001-07-01"), 
             y = 170,
             label = "First period \nof the study",
             hjust = 0.5) +
    geom_vline(xintercept = as.Date("2005-01-01"))+
    geom_vline(xintercept = as.Date("2009-01-01"))+
    annotate(geom = "text", 
             colour = "black", 
             x = as.Date("2007-01-01"), 
             y = 170,
             label = "Second period \nof the study",
             hjust = 0.5) +
    scale_colour_gradient2("MHW cumulative \nintensity", 
                           low = "#ffc866",
                           mid = "#ff6900",
                           high = "#9e0000",
                           limits = c(0, 170))+
    labs(x = NULL, 
         y = expression(atop("Cumulative intensity", 
                             paste("(", degree, "C" %*%"days)"))))+
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,190))+
    xlim(as.Date("1999-01-01"), as.Date("2010-01-01"))+
    theme(legend.title = element_text(hjust = .5))) 

# Join them 

(figureS3 <- p1/p2+plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "a")&
    theme(plot.tag = element_text(face = "bold", size = 14)))

# Save the plot

ggsave("Results/Figure S3.pdf", figureS3,
       width = 10, height = 6)

# Figure S4: Population comparisons --------------------------------------------

load("Data/ResTemp.RData")

(g1 <- dem_cat %>% 
    ggplot(aes(x=cat, y=resistance, fill=cat)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 4,
               position = position_jitter(seed = 2, width = .2),
               shape=21, alpha=0.7) + 
    scale_fill_manual(values = c("#D9AA1E", "#A6122D"))+
    coord_cartesian(xlim = c(1.2, NA),
                    ylim=c(0, 1), 
                    clip = "off") +
    facet_wrap(~site)+
    labs(x="", y="Resistance")+
    theme(legend.position = "none"))

(g2 <- dem_cat %>% 
    ggplot(aes(x=cat, y=speed_recovery, fill=cat)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 4,
               position = position_jitter(seed = 2, width = .2),
               shape=21, alpha=0.7) + 
    facet_wrap(~site)+
    coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    scale_fill_manual(values = c("#D9AA1E", "#A6122D"))+
    labs(x="", y="Speed of recovery")+
    theme(legend.position = "none"))

(figS4 <- g1+g2+plot_annotation(tag_levels = "a"))

ggsave(figS4, 
       filename = "Figure S4.pdf", 
       path = ResultPath,
       width = 8, height = 4)
