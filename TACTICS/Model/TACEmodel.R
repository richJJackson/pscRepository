### TACE Model

## Set working directory
setwd("/Volumes/RICHJ23/Accademic/Synthetic Control/TACTICS")


### Loading functions
source("/Users/richardjackson/Dropbox/Documents/R Utils/statsTools.R")
library(survival)
library(flexsurv)


### Getting Data
tace <- read.csv("pre-TACE model data.csv")

## defining censoring indicator
tace$cens <- as.numeric(tace$death=="Yes")

## defining survival function
tace$s.ob <- Surv(tace$survival..months.+0.1,tace$cens)


### Data Cleaning
tace$tumour.number[which(tace$tumour.number=="")] <- NA
tace$tumour.number <- relevel(factor(tace$tumour.number),ref="Solitary")
table(tace$death)

tace$tum.siz <- log(tace$baseline.maximum.tumour.size..cm.+.05,10)

tace$afp <- log(tace$baseline.AFP..ng.ml.+1,10)

tace$alb <- tace$baseline.albumin..g.l.

tace$bil <- log(tace$baseline.bilirubin..µmol.l.,10)

tace$vi <- tace$vascular.invasion
tace$vi[which(tace$vi=="")] <- NA

tace$aet <- tace$aetiology
tace$aet[which(tace$aet=="")] <- NA
tace$aet <- relevel(factor(tace$aet),ref="HCV")

tace$ecog <- factor(cut(tace$ECOG,c(-0.5,1.5,4),c(0,1)))

tace$cp <- factor(tace$Child.Pugh.grade)
tace$cp[which(tace$cp=="")] <- NA
tace$cp <- as.character(tace$cp)

table(tace$aet)

tace$hcv <- 0
tace$hcv[tace$aet=="HCV"]<-1

tace$hbv <- 0
tace$hbv[tace$aet=="HBV"]<-1

table(tace$hcv,tace$hbv)
### Selecting training set
train <- tace[which(tace$groups.after.random.split..1.training.set..2.internal.validation.set.==1),]
train[1:3,]
train$vascular.invasion

boxplot(train$bil~train$vi)

cm.mod <- coxph(s.ob~tumour.number+tum.siz+afp+vi+alb+bil+hcv+hbv+ecog+cp,data=train)
summary(cm.mod)
anova(cm.mod)

fpm.tace <- flexsurvspline(s.ob~tumour.number+tum.siz+afp+alb+bil+hcv+hbv+vi+ecog+cp,k=1,data=train)

fpm.tace


save(fpm.tace,file="fpm.tace.R")
plot(fpm.tace)




