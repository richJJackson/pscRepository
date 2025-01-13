### PSC Analysis

#Saving 'Home'
hm <- getwd()

## Setting working directory for data
setwd("Data")

### Loading packaged
library(psc)

### Loading data
load("analData.R")


### Counter Factual Model
data_or  <- data[which(data$t2=="OR"),] # Selecting only open repair patients
or_model <- glm(out~I_sex+diameter+asa_cat+(age)+htn+neck_lgth+bmi+i_crea,data=data_or,family="binomial")

summary(or_model)
anova(model)


### Comparing EVAR to OR using psc
data_evar <- data[which(data$t2=="EV"),]

psc_ev <- pscfit(or_model,data_evar)


summary(psc_ev)
plot(psc_ev)



