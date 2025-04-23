
###### Code detialing the process for creating a CFM object from a model,
### UNdertaking model validation and saving structures for sharing on github


### L
library(knitr)
library(RColorBrewer)
library(ggplot2)
library(waffle)
library(ggpubr)
library(gtsummary)



### load psc
setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()


### Loading model
setwd("/Volumes/richj23/Projects/Cancer/HCC/EURAB/May25")

load("AtezoBevmodel_with_types.R")



### Creating model
CFM <- pscCFM(flexspline_model1)
save(CFM,file="atezoBev_CFM.Rds")

CFM$datasumm


### Creating CFM object




