# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        20/07/2023
# - DESCRIPTION: Code create the figures of the manuscript
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #
# Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 
CodePath <- paste0(path,"/Code")

# Libraries

library(dplyr)
library(tidyverse)
library(patchwork)

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

# Figure 2: SNA-LTRE Analysis --------------------------------------------------

load(paste0(ResultPath, "/LTRE.RData"))

# Estimate the porportion of the contribution 

ltre_dat_abs %>% 
  mutate(total=sum(abs(cont))) %>% 
  group_by(component) %>% 
  mutate(cont=sum(abs(cont))/total)

# First we plot the absolute contributions to lambda

(g1 <- ltre_dat_abs %>% 
   ggplot(aes(x=component, y=cont, fill = cont >= 0))+
   geom_bar(stat = "identity")+
   geom_hline(yintercept = 0)+
   scale_fill_manual(values = c("grey70", "grey30"))+
   labs(y=expression(paste("Contributions to ",Delta, " ", lambda[s])), 
        x="")+
   ylim(-0.12, 0.12)+
   scale_x_discrete(labels=c(expression(paste(Delta, " correlation")),
                             expression(paste(Delta, " CV")),
                             expression(paste(Delta," elasticity")),
                             expression(paste(Delta, " means"))))+
   theme(legend.position = "none"))

# Contributions of the demographic processes

(g2 <- ltre_dat_mean %>% 
    mutate(demogr=gsub("Growth", "Progression", demogr)) %>% 
    ggplot(aes(x=demogr, y=cont_mean, fill = cont_mean >= 0))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c("grey70", "grey30"))+
    geom_hline(yintercept = 0)+
    labs(y=expression(atop(paste("Contributions to ",Delta, " ", lambda[s]), 
                           paste("by ", Delta, "means"))), 
         x="Demographic process")+
    ylim(-0.06, 0.06)+
    theme(axis.text.x = element_text(angle = 25, 
                                     hjust = 1, 
                                     vjust=1),
          legend.position = "none"))

# Combine the figures 

(fig2 <- (g1+g2)+plot_annotation(tag_levels = "a"))

# Save the plot 

ggsave("Figure 2.pdf", fig2,
       path = ResultPath,
       width = 12, height = 4)


# Figure 3: Resilience MHWs vs nMHWs -------------------------------------------

load(paste0(DataPath, "/ResTemp.RData"))

# Resistance 

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
    labs(x="", y="Resistance")+
    theme(legend.position = "none"))

# Time of recovery

(g2 <- dem_cat %>% 
    ggplot(aes(x=cat, y=speed_recovery, fill=cat)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 4,
               position = position_jitter(seed = 5, width = .2),
               shape=21, alpha=0.7) + 
    coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    scale_fill_manual(values = c("#D9AA1E", "#A6122D"))+
    labs(x="", y="Speed of recovery")+
    theme(legend.position = "none"))

# Rearrange the dataframe so we just focus on the elasticities

sens_data <- dem_cat %>% ungroup %>% 
  dplyr::select(cat, e_res_fecundity,e_res_growth,e_res_shrinkage,e_res_stasis,
         e_rec_fecundity, e_rec_growth, e_rec_shrinkage, e_rec_stasis) %>% 
  pivot_longer(cols=2:9, names_to = "elasticity", values_to = "value") %>% 
  mutate(resilience=case_when(
    grepl("_res_", elasticity) ~ "resistance",
    grepl("_rec_", elasticity) ~ "recovery",
    TRUE ~ "other"),
    vital_rate=case_when(
      grepl("_fecundity", elasticity) ~ "Fecundity",
      grepl("_growth", elasticity) ~ "Progression",
      grepl("_shrinkage", elasticity) ~ "Shrinkage",
      grepl("_stasis", elasticity) ~ "Stasis",
      TRUE ~ "other"))

# Sensitivity to resistance 

(g3 <- sens_data %>% 
    filter(resilience=="resistance") %>% 
    ggplot(aes(x=vital_rate, y=value, 
               fill=cat)) +
    geom_hline(yintercept = 0, lty="dashed", color="grey25")+
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group=cat),
               size = 2,
               position = position_jitterdodge(jitter.width = .2),
               shape=21, alpha=0.7) + 
    scale_fill_manual("", values = c("#D9AA1E", "#A6122D"))+
    labs(y="Resistance elasticities", 
         x="Demographic process")+
    theme(legend.position = "none"))

# Elasticity to speed of recovery

(g4 <- sens_data %>% 
    filter(resilience=="recovery") %>% 
    ggplot(aes(x=vital_rate, y=value, 
               fill=cat)) +
    geom_hline(yintercept = 0, lty="dashed", color="grey25")+
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(group=cat),
               size = 2,
               position = position_jitterdodge(jitter.width = .2),
               shape=21, alpha=0.7) + 
    scale_fill_manual("", values = c("#D9AA1E", "#A6122D"))+
    labs(y="Speed of recovery \nelasticities", x="Demographic process")+
    theme(legend.position = "none"))

# Join the plots

(top <- g1+g2)

(bottom <- g3+g4+
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom"))

(fig3 <- top/bottom+plot_annotation(tag_levels = "a")&
    theme(legend.text = element_text(size = 14)))

# Save it

ggsave("Figure 3.pdf", fig3, 
       width = 11, height = 8, path = ResultPath)

# Figure 4: Projections --------------------------------------------------------

load(paste0(ResultPath, "/Simulations.RData"))

# Parameters of the simulation

probs <- seq(0, 1, by=0.2)# Probabilities of disturbance
tmax <- 100 # Time of projection
reps <- 1000 # Number of simulation

# Create a gradient of colours 

pal <-  colorRampPalette(c("#D9AA1E", "#A6122D"))
cols <- pal(length(probs))

# Find the minimum time when population goes extinct

(position <- results %>% 
    mutate(probs=paste(prob*100, "MHWs over 100 years"),
           probs=gsub("100 MHWs over 100 years", "Yearly", probs),
           probs=gsub("\\b0 MHWs over 100 years\\b", 
                      "Never", probs),
           probs=fct_reorder(probs, prob)) %>%
    group_by(probs, rep, time) %>% 
    mutate(above.qe= n >= 10,
           alife=cumprod(above.qe)) %>% 
    group_by(probs, time) %>% 
    summarise(prob= 1-(sum(alife, na.rm = T)/reps)) %>% 
    filter(prob==1) %>% 
    group_by(probs) %>% 
    summarise(time=min(time)))


# Quasi-extinction 

(p1<- results %>%
    mutate(probs=paste(prob*100, "MHWs over 100 years"),
           probs=gsub("100 MHWs over 100 years", "Yearly", probs),
           probs=gsub("\\b0 MHWs over 100 years\\b", 
                      "Never", probs),
           probs=fct_reorder(probs, prob)) %>%
    group_by(probs, rep, time) %>% 
    mutate(above.qe= n >= 10,
           alife=cumprod(above.qe)) %>% 
    group_by(probs, time) %>% 
    summarise(prob= 1-(sum(alife, na.rm = T)/reps)) %>% 
    ggplot(aes(x=time, y=prob, colour=probs))+
    geom_line() +
    scale_colour_manual("MHWs frequency",
                        values=cols)+
    labs(y= "Quasi-extinction probability",
         x="Time")+
    xlim(0, 75)+
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
      legend.position = c(.95, .65),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)))


# Resistance 

(p2 <- results %>%
    mutate(probs=as.factor(prob*100),
           probs=gsub("\\b0\\b", "Never", probs),
           probs=gsub("\\b100\\b", "Yearly", probs),
           probs=fct_reorder(probs, prob)) %>% 
    ggplot(aes(x=probs, y=resistance, fill=probs))+
    geom_boxplot(width = .7,
                 outlier.shape = NA,
                 colour="black") +
    theme(legend.position = "none")+
    scale_fill_manual(values=cols)+
    scale_colour_manual(values=cols)+
    ylim(0,1)+
    labs(y= "Resistance",
         x="Frequency (number of MHWs \nover 100 years)"))

# Time of recovery 

(p3 <- results %>%
    mutate(probs=as.factor(prob*100),
           probs=gsub("\\b0\\b", "Never", probs),
           probs=gsub("\\b100\\b", "Yearly", probs),
           probs=fct_reorder(probs, prob)) %>% 
    ggplot(aes(x=probs, y=recovery, fill=probs))+
    geom_boxplot(width = .7,
                 outlier.shape = NA,
                 colour="black") +
    theme(legend.position = "none")+
    scale_fill_manual(values=cols)+
    scale_colour_manual(values=cols)+
    ylim(1,1.2)+
    labs(y= "Speed of recovery",
         x="Frequency (number of MHWs \nover 100 years)"))

# Combine figures 

(fig4 <- p1/(p2+p3)+plot_annotation(tag_levels = "a"))

# Save the figures

ggsave(fig4, 
       filename = "Figure 4.pdf", 
       path = ResultPath,
       width = 10, height = 9)

