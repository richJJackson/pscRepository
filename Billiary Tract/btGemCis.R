

#### Predictive Model for GemCis

### Getting Data
setwd("/Volumes/RICHJ23/Projects/Billiary Tract")
load("modelData.R")

### loading packages
library(flexsurv)

data <- data[-which(data$ps==3),]

datGC <- data[data$Trt=="GemCis",]


### Training and Validation set
set.seed(21032019)
smp <- sample(nrow(datGC),700)
datGC_train <- datGC[smp,]
datGC_valid <- datGC[-smp,]

### Model
fpm <- flexsurvspline(s.ob~diff+ps+primary+livMet+hb_dum/hb2+ace_dum/lace2+alb_dum/alb2+ca_dum/lca1992+ ne_dum/lneut2,data=datGC_train,k=1)


### Validation
fpm_valid <- flexsurvspline(s.ob~diff+ps+primary+livMet+hb_dum/hb2+ace_dum/lace2+alb_dum/alb2+ca_dum/lca1992+ ne_dum/lneut2,data=datGC_valid,k=1)

co <- coef(fpm)
est <- co[-c(1:3)]
gam <- co[1:3]

mm_train <- model.matrix(fpm)
mm_valid <- model.matrix(fpm_valid)

pred_train <- (est)%*%t(mm_train)
pred_valid <- (est)%*%t(mm_valid)

data.fpm <- data.frame(fpm$data[[1]])
data.valid <- data.frame(fpm_valid$data[[1]])


### creating risk groups
p_quant <- quantile(pred_train,c(0.15,0.5,0.85))
rsk_grp_train <- cut(pred_train,c(-Inf,p_quant,Inf),c("risk_grp_1","risk_grp_2","risk_grp_3","risk_grp_4"))

rsk_grp_valid <- cut(pred_valid,c(-Inf,p_quant,Inf),c("risk_grp_1","risk_grp_2","risk_grp_3","risk_grp_4"))



### Calculating Survival Function
modp <- function(x){
	x*(sign(x)+1)/2
}

kn <- fpm$knots
lam <- (max(kn)-kn[2])/(max(kn)-min(kn))
time <- seq(1e-1,72,length=1000)
logt <- log(time)
lam

z2 <- modp(logt-kn[2])^3 - lam*modp(logt-kn[1])^3 - (1-lam)*modp(logt-kn[3])^3
H0 <- exp(gam[1]+gam[2]*logt + gam[3]*z2)
S <- exp(-H0)
plot(fpm)
lines(exp(logt),S,col=3,lwd=4)



### Risk Groups
HR_train <- exp(tapply(pred_train,rsk_grp_train,mean))
HR_valid<- exp(tapply(pred_valid,rsk_grp_valid,mean))


### Set working group for deliverables
setwd("/Volumes/RICHJ23/Projects/Billiary Tract/Deliverables/oct2021")

### Kaplan Meier plots
KMplot(data.fpm$time,data.fpm$status,rsk_grp_train,lwd=5,col=c("forestgreen","royalblue","darkorchid2","darkorange"),xlim=c(0,72),main="Training Set",cex.main=3,ylab="Overall Survival")
lines(exp(logt),S^HR_train[1],col=3,lty=3,lwd=3)
lines(exp(logt),S^HR_train[2],col=4,lty=3,lwd=3)
lines(exp(logt),S^HR_train[3],col=6,lty=3,lwd=3)
lines(exp(logt),S^HR_train[4],col=2,lty=3,lwd=3)
dev.copy2pdf(file="gemCis_train.pdf");dev.off()

KMplot(data.valid$time, data.valid$status,rsk_grp_valid,lwd=5,col=c("forestgreen","royalblue","darkorchid2","darkorange"),xlim=c(0,72),main="Validation Set",cex.main=3,ylab="Overall Survival")
lines(exp(logt),S^HR_valid[1],col=3,lty=3,lwd=3)
lines(exp(logt),S^HR_valid[2],col=4,lty=3,lwd=3)
lines(exp(logt),S^HR_valid[3],col=6,lty=3,lwd=3)
lines(exp(logt),S^HR_valid[4],col=2,lty=3,lwd=3)
dev.copy2pdf(file="gemCis_valid.pdf");dev.off()



