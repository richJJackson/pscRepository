####  Fitting Prosash 3 model

setwd("/Users/Richardjackson/Dropbox/Synthetic Control/PROSASH")

### laoding packages
library(mvtnorm)
library(survival)
library(flexsurv)



### loading data
pros <- read.csv("imputed_data_plus_xb.csv")


#### Functions
modp <- function(x){
	x*(sign(x)+1)/2
}

acc <-  function(old,new){
	ret <- FALSE
	r <- runif(1)
	e <-exp(old-new)	
	if(e>r) ret <- TRUE
	ret
}




### Fit the prosash model 

pros$event <- as.numeric(pros$death)-1
pros$age60 <- pros$age-60
pros$aet <- relevel(pros$aet,ref="HCV")

pros[1:3,]
nrow(pros)
pros[1:3,]

fpm <- flexsurvspline(Surv(os,event)~vi/age60+factor(ecog)+ allmets +logafp+alb+logcreat+logast+aet,data=pros,k=1)
save(fpm,file="flexParaSoraf2.R")

fpm <- flexsurvspline(Surv(os,event)~vi/age60+factor(ecog)+ allmets +logafp+alb+logcreat+logast+aet+lognlr,data=pros,k=1)
save(fpm,file="flexParaSoraf3.R")

fpm <- flexsurvspline(Surv(os,event)~vi/age60+factor(ecog) +logafp+alb+logcreat+logast+aet,data=pros,k=1)
save(fpm,file="flexParaSoraf_nonehs.R")