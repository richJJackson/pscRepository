
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

## loading functions
library(readxl)
library(dplyr)
library(survival)
library(flexsurv)
library(survminer)
library(ggplot2)

#library(devtools)

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


tace_data[1:3,]
tace_data$Age<-round(tace_data$Age)
tace_data$time<-tace_data$time/30.44

tace_data1<-filter(tace_data, treatment=="TACE + Placebo")

tace_data2<-tace_data1[,which(names(tace_data1)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Focality","Hepatitis B","Hepatitis C","Lesion1","Lesion2","Vascular.invasion","Child.Pugh.grade","cen","time","treatment"))]

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



### TACE model
#####  backward selection method
tace_data3[1:3,]


tace_data3 <- tace_data3[tace_data3$time > 0, ]

### AFP may work better on the log scale...
tace_data3$AFP <- log(tace_data3$AFP+1)
### Make sure that Cirrhosis is a factor
tace_data3$Cirrhosis <- as.factor(tace_data3$Cirrhosis)


cox_full<-coxph(Surv(time,cen)~Gender+Age+AFP+Cirrhosis+Bilirubin+Albumin+Lesion1,data=tace_data3)
summary(cox_full)

backward_model<-step(cox_full, direction="backward")
summary(backward_model)


flexspline_model1 <- flexsurvspline(Surv(time, cen) ~ AFP+Cirrhosis+Albumin, data = tace_data3,k=1)  
flexspline_model1


### Visualisation tools to look at the dataset
#install.packages("waffle")
library(waffle)
### Visualisation tools to look at the dataset
#install.packages("RColorBrewer")
library(RColorBrewer)
cfm.vis <- cfmDataVis(flexspline_model1)
ggarrange(cfm.vis[[1]],cfm.vis[[2]],cfm.vis[[3]])


validation1<-cfmValid.flexsurvreg(flexspline_model1)
validation1


setting1<-cfmSett(flexspline_model1)
setting1

flexspline_model1
flexspline_model6

### Setting working directory for saving models
### Please note I've added flexspline_model6 as I think this included the terms from the stepwise procedure
#setwd("~/Documents/GitHub/pscLibrary/TACE_SIRT/Data")
save(flexspline_model1, file="1knotmodel.Rdata")

load("1knotmodel.RData")



###PSC comparison
##Compare SIRT data to TACE model for overall comparison
sirt_data1<-sirt_data[,which(names(sirt_data)%in%c("Patient.TNO","Gender","Age","AFP","Cirrhosis",
"Bilirubin","Albumin","Lesion1","Morphology","PortalVeinInvasion","cen","time","treatment"))]
sirt_data1 <- as.data.frame(sirt_data1)

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

sirt_data1<-na.omit(sirt_data1)


# Align factor levels 
#levels(sirt_data1$Cirrhosis) <- levels(flexspline_model1$data$m$Cirrhosis)

library(waffle)
### Visualisation tools to look at the dataset
#install.packages("RColorBrewer")
library(RColorBrewer)

# Run the comparison function
visComp_sirt <- cfmDataComp(flexspline_model1, sirt_data1)

visComp_sirt <- cfmDataComp(flexspline_model1,sirt_data1)
ggarrange(visComp_sirt[[1]],visComp_sirt[[2]],visComp_sirt[[3]])


psc2<-pscfit(flexspline_model1,sirt_data1)
summary(psc2)
plot(psc2)
summary(psc2)
plot_ite(psc2)

### We can check if the new fit looks reasonable (notice beta=-0.915)
s2 <- surv_fpm(psc2$DC_clean,beta=-0.8801)
s_data <- data.frame("time"=s2$time,"S"=s2$S)
p <- plot(psc2)
p+ geom_line(data=s_data, aes(time,S),col=4,lwd=1.5)


# Construct ggplot manually
p <- ggplot() +
  geom_line(data = s_data, aes(x = time, y = S), col = "pink", lwd = 1.5) +
  labs(x = "Time", y = "Survival Probability") +
  theme_minimal()

print(p)  

### Not a great fit - likely there is some sort of delayed effect (np hazards)

exp(-0.91)



# Summarize SIRT dataset
sirt_summary <- sirt_data1 %>%
  summarise(
    #Age_Mean = mean(Age, na.rm = TRUE),
    #Age_SD = sd(Age, na.rm = TRUE),
    AFP_Mean = mean(AFP, na.rm = TRUE),
    AFP_SD = sd(AFP, na.rm = TRUE),
    Bilirubin_Mean = mean(Bilirubin, na.rm = TRUE),
    Bilirubin_SD = sd(Bilirubin, na.rm = TRUE),
    Albumin_Mean = mean(Albumin, na.rm = TRUE),
    Albumin_SD = sd(Albumin, na.rm = TRUE),
    Lesion1_Mean = mean(Lesion1, na.rm = TRUE),
    Lesion1_SD = sd(Lesion1, na.rm = TRUE),
   #Cirrhosis_Percent = mean(Cirrhosis == "Yes", na.rm = TRUE) * 100,
    #Gender_Male_Percent = mean(Gender == "Male", na.rm = TRUE) * 100,
    #Dead_Percent = mean(cen == 1, na.rm = TRUE) * 100,
    #Mean_Survival_Time = mean(time, na.rm = TRUE)
  ) %>%
  mutate(Treatment = "SIRT")

library(tidyr)
# Compute percentages for all categorical variables
gender_sirt_percent <- sirt_data1 %>%
  group_by(Gender) %>%
  summarise(Percentage = n() / nrow(sirt_data1) * 100) %>%
  pivot_wider(names_from = Gender, values_from = Percentage, names_prefix = "Gender_")

cirrhosis_sirt_percent <- sirt_data1 %>%
  group_by(Cirrhosis) %>%
  summarise(Percentage = n() / nrow(sirt_data1) * 100) %>%
  pivot_wider(names_from = Cirrhosis, values_from = Percentage, names_prefix = "Cirrhosis_")


morphology_percent <- sirt_data1 %>%
  group_by(Morphology) %>%
  summarise(Percentage = n() / nrow(sirt_data1) * 100) %>%
  pivot_wider(names_from = Morphology, values_from = Percentage, names_prefix = "Morphology_")

portalveininvasion_percent <- sirt_data1 %>%
  group_by(PortalVeinInvasion) %>%
  summarise(Percentage = n() / nrow(sirt_data1) * 100) %>%
  pivot_wider(names_from = PortalVeinInvasion, values_from = Percentage, names_prefix = "PortalVeinInvasion_")

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
    Lesion1_Mean = mean(Lesion1, na.rm = TRUE),
    Lesion1_SD = sd(Lesion1, na.rm = TRUE),
    Lesion2_Mean = mean(Lesion2, na.rm = TRUE),
    Lesion2_SD = sd(Lesion2, na.rm = TRUE),
    #Cirrhosis_Percent = mean(Cirrhosis == "Yes", na.rm = TRUE) * 100,
    #Gender_Male_Percent = mean(Gender == "Male", na.rm = TRUE) * 100,
    #Dead_Percent = mean(cen == 1, na.rm = TRUE) * 100,
   # Mean_Survival_Time = mean(time, na.rm = TRUE)
  ) %>%
  mutate(Treatment = "TACE + Placebo")

library(tidyr)

# Compute percentages for all categorical variables
gender_percent <- tace_data3 %>%
  group_by(Gender) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = Gender, values_from = Percentage, names_prefix = "Gender_")

cirrhosis_percent <- tace_data3 %>%
  group_by(Cirrhosis) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = Cirrhosis, values_from = Percentage, names_prefix = "Cirrhosis_")

hepatitisB_percent <- tace_data3 %>%
  group_by(`Hepatitis B`) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = `Hepatitis B`, values_from = Percentage, names_prefix = "HepatitisB_")

hepatitisC_percent <- tace_data3 %>%
  group_by(`Hepatitis C`) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = `Hepatitis C`, values_from = Percentage, names_prefix = "HepatitisC_")

alcohol_percent <- tace_data3 %>%
  group_by(Alcohol) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = Alcohol, values_from = Percentage, names_prefix = "Alcohol_")

focality_percent <- tace_data3 %>%
  group_by(Focality) %>%
  summarise(Percentage = n() / nrow(tace_data3) * 100) %>%
  pivot_wider(names_from = Alcohol, values_from = Percentage, names_prefix = "Focality_")


