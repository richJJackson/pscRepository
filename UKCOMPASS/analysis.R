### Setwd()
setwd("~/Documents/GitHub/pscLibrary/UKCOMPASS")


## loading up packages
library(psc)
library(survival)

## loading data and model
data <- read.csv("Data/data3.csv")
load("Model/5knotmodel 1.Rdata")



### Defining the model

model <- flexspline_model5


### Setting
P <- "Patients undergoing surgery for a complex aortic aneurysm"
I <- "Open Surgery"
C <- "Not Included in Model"
O <- "Overall Survival"
sett <- cfmSett(P,I,C,O)

cfm.ob <- pscCFM(model,setting=sett)


##### PSC

ukc.psc <- pscfit(model,data)


### subgroup for asa = 3,4
id <- which(data$asa%in%c(3,4))
ukc.psc.asa.1 <- pscfit(model,data,id=id)

id <- which(data$asa%in%c(1,2))
ukc.psc.asa.2 <- pscfit(model,data,id=id)

plot(density(ukc.psc$posterior$beta),col=1,lwd=3)
lines(density(ukc.psc.asa.1$posterior$beta),col=2,lwd=3)
lines(density(ukc.psc.asa.2$posterior$beta),col=3,lwd=3)


####

