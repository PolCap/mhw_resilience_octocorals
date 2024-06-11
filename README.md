
# Mediterranean octocoral populations exposed to marine heatwaves are less resilient to disturbances 

Pol Capdevila<sup>1,2</sup>*, Yanis Zentner<sup>1,2</sup>, Graciel·la Rovira<sup>1,2</sup>, Joaquim Garrabou<sup>3,4</sup>, Alba Medrano<sup>1,2</sup>, Cristina Linares<sup>1,2</sup>
 
<sup>1</sup>Departament de Biologia Evolutiva, Ecologia i Ciències Ambientals, Facultat de Biologia, Universitat de Barcelona, Spain. 
<sup>2</sup>Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona (UB), Barcelona, Spain.
<sup>3</sup>Institut de Ciències del Mar-CSIC, Passeig Marítim de la Barceloneta 37-49 08003 Barcelona, Spain.
<sup>4</sup>Aix Marseille Univ, Université de Toulon, CNRS, IRD, MIO, Marseille, France.


#### Contact: pcapdevila[at]ub.edu

---

## Abstract

1. The effects of climate change are now more pervasive than ever. Marine ecosystems have been particularly impacted by climate change, with Marine Heatwaves (MHWs) being a strong driver of mass mortality events. Even in the most optimistic greenhouse gas emission scenarios, MHWs will continue to increase in frequency, intensity, and duration. For this reason, understanding the resilience of marine species to the increase of MHWs is crucial to predicting their viability under future climatic conditions. 

2. In this study, we explored the consequences of MHWs on the resilience (the ability of a population to resist and recover after a disturbance) of a Mediterranean key octocoral species, _Paramuricea clavata_, to further disturbances to their population structure. To quantify _P. clavata_’s capacity to resist and recover from future disturbances, we used demographic information collected from 1999 to 2022, from two different sites in the NW Mediterranean Sea to calculate the transient dynamics of their populations. 

3. Our results showed that the differences in the dynamics of populations exposed and those not exposed to MHWs were driven mostly by differences in mean survivorship and growth. We also showed that after MHWs _P. clavata_ populations had lower resistance and slower rates of recovery than those not exposed to MHWs. Populations exposed to MHWs had lower resistance elasticity to most demographic processes compared to unexposed populations. In contrast, the only demographic process showing some differences when comparing the speed of recovery elasticity values between populations exposed and unexposed to MHWs was stasis. Finally, under scenarios of increasing frequency of MHWs, the extinction of _P. clavata_ populations will accelerate and their capacity to resist and recover after further disturbances will be hampered. 

4. Overall, these findings confirm that future climatic conditions will make octocoral populations even more vulnerable to further disturbances. These results highlight the importance of limiting local impacts on marine ecosystems to dampen the consequences of climate change.


---

## Data

- __`PortCros.csv`__: demographic survey data for Port-Cros National Park populations. 
- __`Medes.csv`__: demographic survey data for Montgrí, Medes Islands and Baix Ter Natural Park populations. 
- __`Matrices.RData`__: list of matrices from Port-Cros National Park and Montgrí, Medes Islands and Baix Ter Natural Park populations. 
- __`DemRes.RData`__: demographic resilience estimates for Port-Cros National Park and Montgrí, Medes Islands and Baix Ter Natural Park populations. 
- __`temperature.RData`__: MHW estimates for the Port-Cros National Park and Montgrí, Medes Islands and Baix Ter Natural Park populations. 

The data files T_Medes.txt and T_Port_Cros.csv appear in the scripts but are not stored in this repository. They are available upon request at https://t-mednet.org/.

---

# Code

To run the statistical analyses we used different R scripts. We have numerated them to show the order that we recommend to follow to produce the results: 

- __`1.Matrices.R`__: code to produce the matrix population models from the raw demographic data.
- __`2.DemoResilience.R`__: code to calculate demographic resilience from the matrix population models.
- __`3.Temperature.R`__: code to estimate the different MHWs estimates from the temperature data.
- __`4.ResilienceTemperature.R`__: code to the temperature and demographic resilience estimates.
- __`5.LTRE.R`__: Code to calculate LTRE of _Paramuricea clavata_ populations exposed and not exposed to MHWs.
- __`6.Projections.R`__: Code to project the populations under different frequencies of MHWs.
- __`7.Figures.R`__: Code to create the figures of the paper.
 - __`8.SupplementaryFigures.R`__: Code to produce the supplementary figures of the manuscript.
 - __`DemoFunctions.R`__: Some functions necessary for the analyses.
  
---

# Software

_R version 4.2.3 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .

Cite the code: [![DOI](https://zenodo.org/badge/436702812.svg)](https://zenodo.org/badge/latestdoi/436702812)
