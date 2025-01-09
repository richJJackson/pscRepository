### Code for data cleaning

## Setting directory
setwd("M:/University/Project/TACE_SIRT1/Data")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.flexsurvreg.R")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.glm.R")
source("M:/University/Project/TACE_SIRT1/functions/plot_ite.R")
source("M:/University/Project/TACE_SIRT1/functions/pscCFM.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmSett.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmValid.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmDataVis.R")
source("M:/University/Project/TACE_SIRT1/functions/cfmValid.flexsurvreg.R")
source("M:/University/R shiny_richard/Richrd_we page/functions.R") 
#source("M:/University/Project/TACE_SIRT1/functions/numVis.R")
#source("M:/University/Project/TACE_SIRT1/functions/numVisComp.R")

## loading functions
library(readxl)
library(dplyr)
library(survival)
library(flexsurv)
library(survminer)
library(ggplot2)

#library(devtools)
#devtools::install_github("richJJackson/psc")

library(psc)
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

tace_data$Age<-round(tace_data$Age)
tace_data$time<-tace_data$time/30.44

tace_data1<-filter(tace_data, treatment=="TACE + Placebo")

tace_data2<-tace_data1[,which(names(tace_data1)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Focality","Lesion1","Vascular.invasion","Child.Pugh.grade","cen","time","treatment"))]

tace_data3<-na.omit(tace_data2)


names(sirt_data)[2]<-"Gender"
names(sirt_data)[6]<-"cen"
names(sirt_data)[9]<-"Albumin"
names(sirt_data)[11]<-"Bilirubin"
names(sirt_data)[14]<-"Lesion1"


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




km3<-survfit(Surv(time,cen)~Gender, data=tace_data3)
km3
summary(km3)

cox3<-coxph(Surv(time,cen)~Gender, data=tace_data3)
summary(cox3)

plot<-ggsurvplot(km3, 
  data = tace_data3,
  conf.int = TRUE,        # Add confidence intervals
  pval = TRUE,            # Display p-value
  risk.table = TRUE,      # Add risk table
  risk.table.y.text.col = FALSE, 
  risk.table.fontsize = 3,
  tables.height = 0.3,
  xlab = "Time (months)", 
  ylab = "Survival probability (%)",
  surv.scale = "percent", 
  risk.table.y.text = TRUE, # Ensure the risk table labels are adjustable
  title = "Kaplan-Meier Survival Curve",
  palette = c("blue", "red"))


    # Customize y-axis percentage display without '%'
    plot$plot <- plot$plot +
      scale_x_continuous(limits = c(0, max(tace_data3$time, na.rm = TRUE)), 
                         expand = c(0, 0.5)) +  # Ensure the x-axis starts from 0
      scale_y_continuous(labels = function(x) as.integer(x * 100)) +
           
      annotate(
        "text", 
        x = 30, y = 0.80, # Position annotation appropriately
        label = "HR = 0.75 (95% CI: 0.39, 1.42)", 
        size = 4, 
        hjust = 0
      )
    
print(plot)




km4<-survfit(Surv(time,cen)~Age, data=tace_data3)
km4

km5<-survfit(Surv(time,cen)~Bilirubin, data=tace_data3)
km5

km6<-survfit(Surv(time,cen)~Albumin, data=tace_data3)
km6


#####  KM plots for overall survival sirt survival


km10 <- survfit(Surv(time, cen) ~ 1, data = sirt_data)
km10

# Plot the Kaplan-Meier curve with ggsurvplot

plot(km10, 
  xlab = "Time (months)", 
  ylab = "Survival Probability", 
  main = "Kaplan-Meier Survival Curve", 
  col = "blue", 
  lwd = 2)

# Add grid for better visualization
grid()


ggsurvplot(km10, 
  data = sirt_data,
  conf.int = TRUE,        # Add confidence intervals
  pval = TRUE,            # Display p-value
  risk.table = TRUE,      # Add risk table
  xlab = "Time (months)", 
  ylab = "Survival Probability", 
  title = "Kaplan-Meier Survival Curve",
  palette = "blue")




km30<-survfit(Surv(time,cen)~Gender, data=sirt_data)
km30
summary(km30)

cox30<-coxph(Surv(time,cen)~Gender, data=sirt_data)
summary(cox30)

plot1<-ggsurvplot(km30, 
  data = sirt_data,
  conf.int = TRUE,        # Add confidence intervals
  pval = TRUE,            # Display p-value
  risk.table = TRUE,      # Add risk table
  risk.table.y.text.col = FALSE, 
  risk.table.fontsize = 3,
  tables.height = 0.3,
  xlab = "Time (months)", 
  ylab = "Survival probability (%)",
  surv.scale = "percent", 
  risk.table.y.text = TRUE, # Ensure the risk table labels are adjustable
  title = "Kaplan-Meier Survival Curve",
  palette = c("blue", "red"))


    # Customize y-axis percentage display without '%'
    plot1$plot <- plot1$plot +
      scale_x_continuous(limits = c(0, max(sirt_data$time, na.rm = TRUE)), 
                         expand = c(0, 0.5)) +  # Ensure the x-axis starts from 0
      scale_y_continuous(labels = function(x) as.integer(x * 100)) +
           
      annotate(
        "text", 
        x = 30, y = 0.80, # Position annotation appropriately
        label = "HR = 0.81 (95% CI: 0.52, 1.27)", 
        size = 4, 
        hjust = 0
      )
    
print(plot1)










### TACE model
#####  backward selection method

cox_full<-coxph(Surv(time,cen)~Gender+Age+AFP+Cirrhosis+Bilirubin+Albumin+Focality+Lesion1,data=tace_data3)
summary(cox_full)

backward_model<-step(cox_full, direction="backward")
summary(backward_model)

# Forward selection
#  fit a null model
null_model <- coxph(Surv(time, cen) ~ 1, data = tace_data3) 

# Then forward selection
forward_model <- step(null_model, direction = "forward", scope = formula(cox_full))
summary(forward_model)

tace_data3 <- tace_data3[tace_data3$time > 0, ]



flexspline_model1 <- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin+Lesion1, data = tace_data3,k=1)  
flexspline_model1


flexspline_model11 <- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin+Focality+Lesion1, data = tace_data3,k=1)  
flexspline_model11

flexspline_model2<- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin+Lesion1, data = tace_data3,k=2)  
flexspline_model2

flexspline_model21 <- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin+Focality+Lesion1, data = tace_data3,k=2)  
flexspline_model21


flexspline_model5 <- flexsurvspline(Surv(stime, dead) ~ AFP+Cirrhosis+Albumin+Lesion1, data = tace_data3,k=5)  
flexspline_model5


flexspline_model6 <- flexsurvspline(Surv(stime, dead) ~ AFP+Cirrhosis+Albumin+Focality, data = tace_data3,k=6)  
flexspline_model6

validation1<-cfmValid.flexsurvreg(flexspline_model1)
validation1

setting1<-cfmSett(flexspline_model1)
setting1



validation2<-cfmValid.flexsurvreg(cfm=flexspline_model1,exData=sirt_data1)
validation2

validation3<-cfmValid.flexsurvreg(flexspline_model1,exData=tace_data13)
validation3

save(flexspline_model1, file="1knotmodel.Rdata")

load("1knotmodel.RData")


### Create on the TACE data using appropriate methodology (stepwise???) to the TACE 2 data
## Once the model has been created using stepwise regression, it can be applied to the TACE 2 dataset(i guess it is Placebo + Sorafenib).

tace_data11<-filter(tace_data,treatment=="TACE + Sorafenib")

tace_data12<-tace_data11[,which(names(tace_data11)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Focality","Lesion1","Vascular.invasion","Child.Pugh.grade","cen","time","treatment"))]

tace_data13<-na.omit(tace_data12)

##comapre with treatment=="TACE + Sorafenib" within Tace dataset using pscfit
psc_tace2<-pscfit(flexspline_model1,tace_data13)
summary(psc_tace2)
plot(psc_tace2)
plot_ite(psc_tace2)


##comapre with treatment=="TACE + Placebo" within Tace dataset using pscfit (the data used for the model itself)
psc_tace1<-pscfit(flexspline_model1,tace_data1)
summary(psc_tace1)
plot(psc_tace1)
plot_ite(psc_tace1)

#### Validate using bootstrapping



###PSC comparison
##Compare SIRT data to TACE model for overall comparison
sirt_data1<-sirt_data[,which(names(sirt_data)%in%c("CUN_ID","Gender","cen","AFP","Albumin","Bilirubin","Cirrhosis","Lesion1","time","Age"))]
psc2<-pscfit(flexspline_model1,sirt_data1)
summary(psc2)

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