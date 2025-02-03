### Code for data cleaning

## Setting directory
setwd("M:/University/Project/TACE_SIRT1/Data")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.flexsurvreg.R")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.glm.R")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.R")
source("M:/University/Project/TACE_SIRT1/functions/pscCFM.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmSett.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmValid.R")
#source("M:/University/Project/TACE_SIRT1/functions/cfmDataVis.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmValid.flexsurvreg.R")
source("M:/University/R shiny_richard/Richrd_we page/functions.R") 
source("M:/University/Project/TACE_SIRT1/functions/numVis.R")
source("M:/University/Project/TACE_SIRT1/functions/numVisComp.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmDataVis (1).R")
source("M:/University/Project/TACE_SIRT1/functions/cfmDataComp.R")
source("M:/University/Project/TACE_SIRT1/functions/facVis.R")
source("M:/University/Project/TACE_SIRT1/functions/facVisComp.R")
#library(devtools)
#devtools::install_github("richJJackson/psc")
#install.packages("waffle")
library(waffle)
### Visualisation tools to look at the dataset
#install.packages("RColorBrewer")
library(RColorBrewer)
library(psc)


#### ASWATHY - this is just what i use to load the functions (you can ignore this)
setwd("/Users/richardjackson/Documents/GitHub/pscVis/pscVis")
devtools::load_all()
setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()


## loading functions
library(readxl)
library(dplyr)
library(survival)
library(flexsurv)
library(survminer)
library(ggplot2)

## loading data
sirt_data <- read_excel("SIRT_UK_MDTA_10_2024.xlsx")
tace_data <- read.csv("tace_2_data_requested_liverpool_2022_11_24.csv")


##Cleaning datasets
names(tace_data)[2]<-"Gender"
names(tace_data)[3]<-"Age"
names(tace_data)[5]<-"Cirrhosis"
names(tace_data)[6]<-"Hepatitis B"
names(tace_data)[7]<-"Hepatitis C"
names(tace_data)[8]<-"Alcohol"
names(tace_data)[11]<-"Bilirubin"
names(tace_data)[12]<-"Albumin"
names(tace_data)[13]<-"Focality"
names(tace_data)[14]<-"Lesion1"
names(tace_data)[15]<-"Lesion2"
names(tace_data)[18]<-"cen"
names(tace_data)[19]<-"time"
names(tace_data)[20]<-"treatment"

tace_data[1:3,]
tace_data$Age<-round(tace_data$Age)
tace_data$time<-tace_data$time/30.44

#tace_data1<-filter(tace_data, treatment=="TACE + Placebo")

tace_data1 <- tace_data[which(tace_data$treatment=="TACE + Placebo"),]
tace_data2<-tace_data1[,which(names(tace_data1)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Focality","Lesion1","Vascular.invasion","Child.Pugh.grade","cen","time","treatment"))]

tace_data3<-na.omit(tace_data2)

names(sirt_data)[1]<-"Patient.TNO"
names(sirt_data)[2]<-"Gender"
names(sirt_data)[6]<-"cen"
names(sirt_data)[9]<-"Albumin"
names(sirt_data)[11]<-"Bilirubin"
names(sirt_data)[14]<-"Lesion1"
names(sirt_data)[20] <-"treatment"

# Convert "Yes" to 1 and "No" to 0 in the Cirrhosis column
sirt_data$Cirrhosis <- ifelse(sirt_data$Cirrhosis == "yes", 1, 0)
sirt_data$cen <- ifelse(sirt_data$cen == "dead", 1, 0)

sirt_data$DateBirth<-as.Date(sirt_data$DateBirth, format = "%d.%m.%Y")
sirt_data$DateTx<-as.Date(sirt_data$DateTx, format = "%d.%m.%Y")
sirt_data$DateLastVis<-as.Date(sirt_data$DateLastVis, format = "%d.%m.%Y")


## Defining primary outcome (overall survival)
sirt_data$stime_days<-as.numeric(difftime(sirt_data$DateLastVis,sirt_data$DateTx, unit="days"))
sirt_data$time<-sirt_data$stime_days/30.44

sirt_data$age_days<-as.numeric(difftime(sirt_data$DateLastVis,sirt_data$DateBirth, unit="days"))
sirt_data$Age<-round(sirt_data$age_days/365.25)
sirt_data <- sirt_data %>%
  mutate(Treatment = "SIRT")

sirt_data<-na.omit(sirt_data)

####1.Summary table summsing data between TACE and SIRT data cohort


# Summarize SIRT dataset
sirt_summary <- sirt_data %>%
  summarise(
    Age_Mean = mean(Age, na.rm = TRUE),
    Age_SD = sd(Age, na.rm = TRUE),
    AFP_Mean = mean(AFP, na.rm = TRUE),
    AFP_SD = sd(AFP, na.rm = TRUE),
    Bilirubin_Mean = mean(Bilirubin, na.rm = TRUE),
    Bilirubin_SD = sd(Bilirubin, na.rm = TRUE),
    Albumin_Mean = mean(Albumin, na.rm = TRUE),
    Albumin_SD = sd(Albumin, na.rm = TRUE),
    Cirrhosis_Percent = mean(Cirrhosis == "Yes", na.rm = TRUE) * 100,
    Gender_Male_Percent = mean(Gender == "Male", na.rm = TRUE) * 100,
    Dead_Percent = mean(cen == 1, na.rm = TRUE) * 100,
    Mean_Survival_Time = mean(time, na.rm = TRUE)
  ) %>%
  mutate(Treatment = "SIRT")

# Summarize TACE dataset
tace_summary <- tace_data3%>%
  summarise(
    Age_Mean = mean(Age, na.rm = TRUE),
    Age_SD = sd(Age, na.rm = TRUE),
    AFP_Mean = mean(AFP, na.rm = TRUE),
    AFP_SD = sd(AFP, na.rm = TRUE),
    Bilirubin_Mean = mean(Bilirubin, na.rm = TRUE),
    Bilirubin_SD = sd(Bilirubin, na.rm = TRUE),
    Albumin_Mean = mean(Albumin, na.rm = TRUE),
    Albumin_SD = sd(Albumin, na.rm = TRUE),
    Cirrhosis_Percent = mean(Cirrhosis == "Yes", na.rm = TRUE) * 100,
    Gender_Male_Percent = mean(Gender == "Male", na.rm = TRUE) * 100,
    Dead_Percent = mean(cen == 1, na.rm = TRUE) * 100,
    Mean_Survival_Time = mean(time, na.rm = TRUE)
  ) %>%
  mutate(Treatment = "TACE + Placebo")

# Combine summaries into a single table
summary_table <- bind_rows(sirt_summary, tace_summary)



##2. Kaplam-Meier plots for tace data
km1 <- survfit(Surv(time, cen) ~ 1, data = tace_data3)
km1

# Plot the Kaplan-Meier curve with ggsurvplot
plot(km1, 
  xlab = "Time (months)", 
  ylab = "Survival Probability", 
  main = "Kaplan-Meier Survival Curve", 
  col = "blue", 
  lwd = 2)

# Add grid for better visualization
grid()


ggsurvplot(km1, 
  data = tace_data3,
  conf.int = TRUE,        # Add confidence intervals
  pval = TRUE,            # Display p-value
  risk.table = TRUE,      # Add risk table
  xlab = "Time (months)", 
  ylab = "Survival Probability", 
  title = "Kaplan-Meier Survival Curve",
  palette = "blue")



#### OMITTED Other KM plots for brevity

tace_data3 <- tace_data3[tace_data3$time > 0, ]

### AFP may work better on the log scale...
tace_data3$AFP <- log(tace_data3$AFP+1)
### Make sure that Cirrhosis is a factor
tace_data3$Cirrhosis <- as.factor(tace_data3$Cirrhosis)

### Visualisation tools to look at the dataset
#install.packages("RColorBrewer")
library(RColorBrewer)
cfm.vis <- cfmDataVis(flexspline_model1)
ggarrange(cfm.vis[[1]],cfm.vis[[2]],cfm.vis[[3]],cfm.vis[[4]])



### TACE model
flexspline_model1 <- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin+Lesion1, data = tace_data3,k=1)  
plot(flexspline_model1)


### Visualisation tools to look at the dataset
cfm.vis <- cfmDataVis(flexspline_model1)
ggarrange(cfm.vis[[1]],cfm.vis[[2]],cfm.vis[[3]],cfm.vis[[4]])

validation1<-cfmValid.flexsurvreg(flexspline_model1)
validation1

setting1<-cfmSett(flexspline_model1)
setting1

### creating model object
cfm.ob <- pscCFM(flexspline_model1,setting=setting1,valid=validation1)
save(cfm.ob,file="TACEcfm.R")

### Setting working directory for saving models




### Create on the TACE data using appropriate methodology (stepwise???) to the TACE 2 data
## Once the model has been created using stepwise regression, it can be applied to the TACE 2 dataset(i guess it is Placebo + Sorafenib).

#tace_data11[1:3,]
### AFP may work better on the log scale...
tace_data$AFP <- log(tace_data$AFP+1)
### Make sure that Cirrhosis is a factor
tace_data$Cirrhosis <- as.factor(tace_data$Cirrhosis)

library(tidyr)
#tace_data11<-filter(tace_data,tace_data$treatment=="TACE + Sorafenib")
tace_data11 <- tace_data[which(tace_data$treatment=="TACE + Sorafenib"),]

tace_data12<-tace_data11[,which(names(tace_data11)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Focality","Lesion1","Vascular.invasion","Child.Pugh.grade","cen","time","treatment"))]

tace_data13<-na.omit(tace_data12)

### need to take logs in AFP again
#tace_data13$AFP <- log(tace_data13$AFP+1)

##comapre with treatment=="TACE + Sorafenib" within Tace dataset using pscfit

### Before fitting a pscfit model we can make sure the data match...
#tace_data13$Cirrhosis<-as.factor(tace_data13$Cirrhosis)


visComp <- cfmDataComp(flexspline_model1,tace_data13)

flexspline_model1$data$m$Albumin <- as.numeric(flexspline_model1$data$m$Albumin)

visComp <- cfmDataComp(flexspline_model1,tace_data13)
ggarrange(visComp[[1]],visComp[[2]],visComp[[3]],visComp[[4]])


psc_tace2<-pscfit(flexspline_model1,tace_data13)
summary(psc_tace2)
p <- plot(psc_tace2)
plot_ite(psc_tace2)

### We can check if the new fit looks reasonable (notice beta=-0.24)
s2 <- surv_fpm(psc_tace2$DC_clean,beta=-0.224)
s_data <- data.frame("time"=s2$time,"S"=s2$S)
p+ geom_line(data=s_data, aes(time,S),col=4,lwd=1.5)


###########################################################################













###PSC comparison
##Compare SIRT data to TACE model for overall comparison


sirt_data1<-sirt_data[,which(names(sirt_data)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Lesion1","cen","time","treatment"))]
#sirt_data1<-sirt_data
#names(sirt_data1)[1]<-"Patient.TNO"
sirt_data1$Patient.TNO<-as.integer(sirt_data1$Patient.TNO)
#### I think the lesion sizes are measured on different scales...
sirt_data1$Lesion1 <- sirt_data1$Lesion1/10
sirt_data1$Cirrhosis <- factor(sirt_data1$Cirrhosis)
sirt_data1$AFP <- log(sirt_data1$AFP+1)

sirt_data1$Gender<-as.character(sirt_data1$Gender)
sirt_data1$cen<-as.integer(sirt_data1$cen)
sirt_data1$Albumin<-as.integer(sirt_data1$Albumin)


sirt_data1$time<-as.numeric(sirt_data1$time)


# Align factor levels 
#levels(sirt_data1$Cirrhosis) <- levels(flexspline_model1$data$m$Cirrhosis)



sirt_data1 <- as.data.frame(sirt_data1)
pscfit(flexspline_model1, sirt_data1)

# Run the comparison function
visComp_sirt <- cfmDataComp(flexspline_model1, sirt_data1)
ggarrange(visComp_sirt[[1]],visComp_sirt[[2]],visComp_sirt[[3]],visComp_sirt[[4]])


psc2<-pscfit(flexspline_model1,sirt_data1)
summary(psc2)
plot(psc2)
summary(psc2)

### We can check if the new fit looks reasonable (notice beta=-0.915)
s2 <- surv_fpm(psc2$DC_clean,beta=-0.915)
s_data <- data.frame("time"=s2$time,"S"=s2$S)
p <- plot(psc2)
p+ geom_line(data=s_data, aes(time,S),col=4,lwd=1.5)

### Not a great fit - likely there is some sort of delayed effect (np hazards)

exp(-0.91)

### Try some sub-group analysis

sirt_data[1:3,]
id <- which(sirt_data1$Albumin<37.5)
psc2_sg1<-pscfit(flexspline_model1,sirt_data1,id=id)
summary(psc2_sg1)
plot(psc2_sg1)
plot_ite(psc2_sg1)


id <- which(sirt_data1$Albumin>37.5)
psc2_sg2<-pscfit(flexspline_model1,sirt_data1,id=id)
summary(psc2_sg2)
plot(psc2_sg2)




p <- cfmVis[[i]]
x <- x

numVisComp <- function(p,x){
  dnew <- data.frame("xnew"=x);dnew
  names(dnew) = "xnew"
  p + geom_density(aes(x=xnew,y=..density..),data=dnew, fill="#404080",color="#404080" )
}

dnew

ggarrange(cfm.vis[[1]],cfm.vis[[2]],cfm.vis[[3]],cfm.vis[[4]])


library(RColorBrewer)


flexspline_model1
sirt_data1[1:3,]
tace_data1[1:3,]

flexspline_model1
str(psc2)
psc2$model.type

plot(psc2)
plot_ite(psc2)

#plot_ite(flexspline_model1)


#pscCFM(cfm=flexspline_model1,covnm=NULL,setting=NULL,valid=NULL,citation=NULL,plot=F)
#numVis(psc2)

#pscCFM(psc2)

psc3<-pscfit(flexspline_model2,sirt_data1)
summary(psc3)
plot(psc3)
plot_ite(psc3)

psc4<-pscfit(flexspline_model2,sirt_data1)
summary(psc4)
plot(psc4)




#### Example data and model from psc package

bin.mod <- psc::bin.mod
data <- psc::data
bin.psc <- pscfit(bin.mod,data)
 plot_ite(bin.psc)


 str(bin.psc)
 bin.psc$model.type