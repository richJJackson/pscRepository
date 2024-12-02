### Synthetic Controls - ESPAC3 and ESPAC4


### laoding packages
library(mvtnorm)
library(survival)
library(flexsurv)
library(knitr)

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

### Getting Data
setwd("/Users/richardjackson/Dropbox/Synthetic Control/ESPAC")

espac3 <- read.csv("espac3.csv")
espac4 <- read.csv("espac4.csv")



########################
########################
########### Fitting ESPAC3 Model

### Calculating OS from surgery date
op_dt <- as.Date(espac3$dop,"%d%b%Y")
ce_dt <- as.Date(espac3$cens_dt,"%d%b%Y")
espac3$stime <- as.numeric(ce_dt-op_dt+1)/30.44
espac3$s.ob <- Surv(espac3$stime,espac3$death)

### Selecting only Gem patients
espac3 <- espac3[which(espac3$trt==1),];nrow(espac3)
espac3 <- espac3[-which(espac3$grade==4),];nrow(espac3)
espac3$diabetes[which(espac3$diabetes==3)] <- 2

### Univariate analysis
espac3$lca199 <- log(as.numeric(as.character(espac3$ca199_post))+1)
espac3$lca199p <- log(as.numeric(as.character(espac3$ca199_pre))+1)

e3 <- espac3[,which(names(espac3)%in%c("s.ob","sex","rstat","t","nodes","m","stage","ps","grade","smoking","diabetes","local_inv","max_size","lca199","lca199p"))]


table(espac3$country)
table(espac3$nodes,espac3$n)
### removing missing data (complete case only for mv build)

miss.id <- which(is.na(e3),arr.ind=T)
miss.id <- unique(miss.id[,1])

e3 <- e3[-miss.id,]

e3$rstat <- as.factor(e3$rstat)
e3$t <- as.factor(e3$t)
e3$nodes <- as.factor(e3$nodes)
e3$m <- as.factor(e3$m)
e3$stage <- as.factor(e3$stage)
e3$ps <- as.factor(e3$ps)
e3$grade <- as.factor(e3$grade)
e3$smoking <- as.factor(e3$smoking)
e3$diabetes <- as.factor(e3$diabetes)
e3$local_inv <- as.factor(e3$local_inv)
e3$max_size <- sqrt(as.numeric(as.character(e3$max_size)))

### removing t=5
e3 <- e3[-which(e3$t==5),]
#e3 <- e3[-which(e3$lca199==0),]
table(e3$t)

### multivariabel
null.mod <- coxph(s.ob~.,data=e3)
step.mod <- step(null.mod,k=3,direction="both")


### Fitting final model
espac3$lca199 <- log(as.numeric(as.character(espac3$ca199_post))+1)


espac3$ca199_post
espac3$lca199
fin.mod <- coxph(s.ob~factor(t)+factor(nodes)+factor(grade)+(lca199),data=espac3)
summary(fin.mod)
anova(fin.mod)

## Creating clean espac3 dataset
espac3 <- espac3[,which(names(espac3)%in%c("s.ob","rstat","nodes","grade","local_inv","t","lca199"))]
espac3 <- espac3[-unique(which(is.na(espac3),arr.ind=T)[,1]),]

### removing t=5 and ca199 = 0 values
espac3 <- espac3[-which(espac3$t==5),]
espac3 <- espac3[-which(espac3 $lca199==0),]


fpm <- flexsurvspline(s.ob~factor(t)+factor(nodes)+factor(grade)+(lca199),data=espac3,k=1)
fpm.res <- fpm$res
fpm.eres <- exp(fpm.res)




r1 <- paste(round(fpm.res[,1],2)," (",round(fpm.res[,4],2),")",sep="")

r2 <- paste(round(fpm.eres[,1],2)," (",round(fpm.eres[,2],2)," - ",round(fpm.eres[,3],2),")",sep="")

fpm.kab <- cbind(r1,r2)
colnames(fpm.kab) <- c("Est (se)","HR (95% CI)")
rownames(fpm.kab) <- rownames(fpm.res)
kable(fpm.kab)

summary(summary(fpm))

save(fpm,file="/Users/richardjackson/Dropbox/Synthetic Control/Models/flexParaGem.R")
save(espac3,file="/Users/richardjackson/Dropbox/Synthetic Control/Models/espac3fpm.R")



### creating covariate dataset and calculating linear predictor
espac3.mm <-model.matrix(fpm)
est <- coef(fpm)
espac3$lp <- espac3.mm%*%est[-c(1:3)]


###
cutP <- c(-Inf,quantile(espac3$lp,c(0.15,0.5,0.85)),Inf)
espac3$riskGrp <- cut(espac3$lp,cutP)


################################
################################
### plotting to get an idea of discrimination in the linear predictor


#### Plotting FPM model
kn <- fpm$knots
lam <- (max(kn)-kn[2])/(max(kn)-min(kn))
gam <- est[1:3]

### Estimation of Hazard function

### building the hazard function from PROSASH2
time <- seq(1e-1,120,length=1000)
logt <- log(time)

### Defining Hazard function
z2 <- modp(logt-kn[2])^3 - lam*modp(logt-kn[1])^3 - (1-lam)*modp(logt-kn[3])^3
H <- exp(gam[1]+gam[2]*logt + gam[3]*z2)
S0 <- exp(-H)

##### basic survival function
S_niave <- S0
risk_haz <- tapply(espac3$lp,espac3$riskGrp,mean)
espac3[1:3,]


sfit <- survfit(s.ob~riskGrp,data=espac3)
plot(sfit,col=c(1,2,3,4,5),lwd=5,xlim=c(0,60),ylab="Overall Survival",xlab="Time (Months)",bty="l")
for(i in 1:5){
	S <- S0^exp(risk_haz[i])
	lines(time,S,col=i,lwd=4,lty=2)	
}

dev.copy2pdf(file="ESPAC3_KM.pdf")
################################
################################









################################
################################
### Validating against ESPAC-4 data

### calculating survival time from oepration
espac4$stime <- as.numeric(as.Date(espac4$cens_date2,"%d%b%Y") - as.Date(espac4$rand_date_chk,"%m/%d/%Y"))/30.44
espac4$s.ob <- Surv(espac4$stime,espac4$event)

## Cleaning ESPAC4 data
espac4$nodes <- cut(espac4$AJCC8Nstage,c(-0.5,1.5,2.5),labels=c(1,2))
espac4$grade <- cut(espac4$differentiation,c(-.5,.5,2.5,4.5),labels=c(1,2,3))
espac4$lca199 <- log(as.numeric(as.character(espac4$ca199_post))+1)

espac4$t <- as.character(espac4$AJCC8Tstage)
espac4$t[which(espac4$t%in%c("T1a","T1b","T1c"))] <- 2
espac4$t[which(espac4$t%in%c("T2"))] <- 3
espac4$t[which(espac4$t%in%c("T3"))] <- 4

espac4$t <- factor(espac4$t,levels=c("2","3","4"))

table(e3$t)
table(espac4$t)
espac4$local_inv <- espac4$local_inv+1
espac4$R_status <- espac4$R_status+1

espac4$lca199 <- log(espac4$ca199_post+1)
espac4$lca199p <- log(espac4$ca199_pre+1)


e4 <- espac4[,which(names(espac4)%in%c("s.ob","nodes","grade","t","R_status","local_inv","lca199","lca199p"))]

names(e4) <- c("grade","local_inv","rstat","s.ob","nodes","lca199","t")

e4gem <- e4[which(espac4$rand_arm_chk==1),]
e4gcap <- e4[which(espac4$rand_arm_chk==2),]

### removing missing data
e4gem <- e4gem[-unique(which(is.na(e4gem),arr.ind=T)[,1]),]
e4gcap <- e4gcap[-unique(which(is.na(e4gcap),arr.ind=T)[,1]),]



#### calculating linear predictor
e4gem.fpm <- flexsurvspline(s.ob~(factor(nodes)+factor(t))+factor(grade)+(lca199),data=e4gem,k=1)
e4gcap.fpm <- flexsurvspline(s.ob~(factor(nodes)+factor(t))+factor(grade)+(lca199),data=e4gcap,k=1)


cbind(coef(fpm),coef(e4gem.fpm),coef(e4gcap.fpm))

### Getting model matrix and calculating linear predictor
mm.e4gem <- model.matrix(e4gem.fpm)
mm.e4gcap <- model.matrix(e4gcap.fpm)

e4gem$lp <- mm.e4gem%*%coef(fpm)[-c(1:3)]
e4gcap$lp <- mm.e4gcap%*%coef(fpm)[-c(1:3)]


################################


stab.cov <- e3[,which(names(e3)%in%c("grade","nodes","lca199","t"))]

espac3$nodes <- as.factor(espac3$nodes)
espac3$t <- as.factor(espac3$t)
espac3$grade <- as.factor(espac3$grade)





stab1 <- summaryTable(espac3[,which(names(espac3)%in%c("grade","nodes","lca199","t"))],rep("Gem (ESPAC-3)",nrow(espac3)),row=FALSE)
stab2 <- summaryTable(e4gem[,which(names(e4gem)%in%c("grade","nodes","lca199","t"))],rep("Gem (ESPAC-4)",nrow(e4gem)),row=FALSE)
stab3 <- summaryTable(e4gcap[,which(names(e4gcap)%in%c("grade","nodes","lca199","t"))],rep("GemCap (ESPAC-4)",nrow(e4gcap)),row=FALSE)

summ.table <- merge(merge(stab1[,-4],stab2[,-4]),stab3[,-4])
write.csv(summ.table,"summTable.csv")






################################
###############
#### Validating ESPAC-4 Gem

##### basic survival function
S_niave <- S0
risk_haz <- tapply(espac3$lp,espac3$riskGrp,mean)

sfit <- survfit(s.ob~riskGrp,data=espac3)
plot(sfit,col=c(1,2,3,4,5),lwd=5,xlim=c(0,60),ylab="Overall Survival",xlab="Time (Months)",bty="l")
for(i in 1:5){
	S <- S0^exp(risk_haz[i])
	lines(time,S,col=i,lwd=4,lty=2)	
}

dev.copy2pdf(file="ESPAC3_KM.pdf")

### Validating using ESPAC 4



### Cox model and median survival
par(mfrow=c(2,1))
hist(espac3$lp,col=4,breaks=20,xlab="",xlim=c(0,4),main="Derivation Data")
abline(v=quantile(espac3$lp,p=c(0.16,0.5,0.84)),lwd=3,col="hotpink",lty=3)

hist(e4gem$lp,col=2,breaks=20,xlab="Prognostic Index (PI)",xlim=c(0,4),main="Validation Data")
abline(v=quantile(e4gem$lp,p=c(0.16,0.5,0.84)),lwd=3,col=5,lty=3)
dev.copy2pdf(file="valid_e4_hist.pdf")




### Calibration Slope
valid.cm <- coxph(s.ob~lp,data=e4gem)
valid.slope <- coef(valid.cm)
valid.slope
valid.con <- round(valid.cm$concordance,3);valid.con

orig.cm <- coxph(s.ob~factor(t)+factor(nodes)+factor(grade)+lca199,data=espac3)


orig.con <- round(orig.cm $concordance,3);orig.con
valid.con <- round(valid.cm$concordance,3)

#        D     se(D)       R.D      R.PM 
#0.8137291 0.1018472 0.1365002 0.1347023 
royston(valid.cm)
#        D     se(D)       R.D      R.PM 
#0.8137291 0.1018472 0.1365002 0.1347023 



setEPS()
postscript("ESPACvalid.eps",height=12,width=8)


par(mfrow=c(3,1))
##### basic survival function
S_niave <- S0
risk_haz <- tapply(espac3$lp,espac3$riskGrp,mean)

sfit <- survfit(s.ob~riskGrp,data=espac3)
plot(sfit,col=c(1,2,3,4,5),lwd=5,xlim=c(0,60),ylab="Overall Survival",xlab="Time (Months)",bty="l")
for(i in 1:5){
	S <- S0^exp(risk_haz[i])
	lines(time,S,col=i,lwd=4,lty=2)	
}
text(5,0.1,c("a)"),cex=2,font=2)
legend(30,0.99,c("KM estimates (ESPAC-3-Gem)","Gem Predictive Model"),lty=c(1,2),cex=1.2,bty="n",lwd=3)




##### basic survival function
cutP2 <- c(-Inf,quantile(e4gem$lp,p=c(0.15,0.5,0.85)),Inf)
e4gem$riskGrp <- cut(e4gem$lp,cutP)

###
survfit(s.ob~riskGrp,data=espac3)
summary(coxph(s.ob~riskGrp,data=espac3))

survfit(s.ob~riskGrp,data=e4gem)
summary(coxph(s.ob~riskGrp,data=e4gem))
###


risk_haz2 <- exp(c(0,coef(coxph(s.ob~riskGrp,data=espac3))))
risk_haz


sfit <- survfit(s.ob~riskGrp,data= e4gem);sfit
plot(sfit,col=c(1,2,3,4,5),lwd=5,xlim=c(0,60),xlab="Time (Months)",ylab="Overall Survival")

risk_haz <- tapply(e4gem $lp, e4gem $riskGrp,mean)

for(i in 1:5){
	S <- S0^exp(risk_haz[i])
	lines(time,S,col=i,lwd=4,lty=3)	
}
text(5,0.1,c("b)"),cex=2,font=2)

legend(30,0.99,c("KM estimates (ESPAC-4-Gem)","Gem Predictive Model"),lty=c(1,2),cex=1.2,bty="n",lwd=3)



mnHR <- mean(e4gcap$lp)

sfit <- survfit(s.ob~1,data=e4gcap);sfit
plot(sfit,col=c(1,2,3,4,5),lwd=4,xlim=c(0,60),xlab="Time (Months)",ylab="Overall Survival",conf.int=F)
S <- S0^exp(mnHR)
lines(time,S,col=2,lwd=4,lty=3)	
text(5,0.1,c("c)"),cex=2,font=2)

legend(30,0.99,c("KM estimates (ESPAC-4-GemCap)","Gem Predictive Model"),lty=c(1,2),col=c(1,2),cex=1.2,bty="n",lwd=3)

dev.off()
#dev.copy2pdf(file="ESPAC3_4_KM.pdf")


###

survfit(s.ob~riskGrp,data=espac3)
summary(coxph(s.ob~riskGrp,data=espac3))

survfit(s.ob~riskGrp,data=e4gem)
summary(coxph(s.ob~riskGrp,data=e4gem))
################################
################################




### ADD in Synthetic Control as a means of validation!!!!

###############





















################################

############### Synthetic Control
### Set Working directory
setwd("/Users/richardjackson/Dropbox/Synthetic Control/Synthetic Control Example")
#source("/Users/richardjackson/Dropbox/Documents/R Utils/statsTools.R")



################################
############## SC
################################

### loading functions
source("synthContFunc.R")





## Calibration for gem

id <- 1:nrow(e4gem)
###########################

kn <- fpm$knots
lam <- (max(kn)-kn[2])/(max(kn)-min(kn))
est <- coef(fpm)




rand_trial_cm <- coxph(s.ob~(factor(nodes)+factor(t))+factor(grade)+(lca199)+rand_arm_chk,data=espac4)
rand_trial_cm2 <- coxph(s.ob~rand_arm_chk,data=espac4)
summary(rand_trial_cm2)
summary(rand_trial_cm)


cov <- mm.e4gem[id,]
tim <- as.numeric(e4gem$s.ob[id])[1:nrow(cov)];tim
cen <- as.numeric(e4gem$s.ob[id])[(nrow(cov)+1):(2*nrow(cov))];cen
as.numeric(e4gcap$s.ob)

op <- optim(fn=lik,list(beta=0),time=tim,cen=cen,cov=cov,est=est,lam=lam,kn=kn,method="Brent",lower=-10,upper=10,hessian=T)
op
op$par
exp(op$par)
sqrt(1/op$hessian)*1.96


###  Synthetic Control - carful not to set the sims too high - it takes a while!
res_null <- synthComp(fpm,tim,cen,cov,25000)

res_null[1:3,]

exp(quantile(res_null[-c(1:5000),10],p=c(0.025,0.5,0.975)))
nrow(res_null)

id <- 1:nrow(e4gcap)
cov <- mm.e4gcap[id,]
tim <- as.numeric(e4gcap$s.ob[id])[1:nrow(cov)];tim
cen <- as.numeric(e4gcap$s.ob[id])[(nrow(cov)+1):(2*nrow(cov))];cen


op <- optim(fn=lik,list(beta=0),time=tim,cen=cen,cov=cov,est=est,lam=lam,kn=kn,method="Brent",lower=-5,upper=5,hessian=T)
op$par
exp(op$par)
sqrt(1/op$hessian)*1.96




###  Synthetic Control - carful not to set the sims too high - it takes a while!
res <- synthComp(fpm,tim,cen,cov,100000)

thin <- seq(5000,100000,by=10)
quantile(res[thin,10],p=c(0.5,0.025,0.975))
exp(quantile(res[thin,10],p=c(0.5,0.025,0.975)))



##############################################
##############################################
### Sub-group analysis

ln.id <- which(cov[,1]==1)
res.ln.pos <- synthComp(fpm,tim[ln.id],cen[ln.id],cov[ln.id,],100000)

ln.id <- which(cov[,1]==0)
res.ln.neg <- synthComp(fpm,tim[ln.id],cen[ln.id],cov[ln.id,],100000)

t.id <- which(cov[,2]==0&cov[,3]==0)
res.t.2 <- synthComp(fpm,tim[t.id],cen[t.id],cov[t.id,],100000)

t.id <- which(cov[,2]==1)
res.t.3 <- synthComp(fpm,tim[t.id],cen[t.id],cov[t.id,],100000)

t.id <- which(cov[,3]==1)
res.t.4 <- synthComp(fpm,tim[t.id],cen[t.id],cov[t.id,],100000)



g.id <- which(cov[,4]==0&cov[,5]==0)
res.g.1 <- synthComp(fpm,tim[g.id],cen[g.id],cov[g.id,],100000)

g.id <- which(cov[,4]==1)
res.g.2 <- synthComp(fpm,tim[g.id],cen[g.id],cov[g.id,],100000)

g.id <- which(cov[,5]==1)
res.g.3 <- synthComp(fpm,tim[g.id],cen[g.id],cov[g.id,],100000)


c.id <- which(as.numeric(cov[,6])<log(20))
res.c.1 <- synthComp(fpm,tim[c.id],cen[c.id],cov[c.id,],100000)

c.id <- which(as.numeric(cov[,6])>=log(20))
res.c.2 <- synthComp(fpm,tim[c.id],cen[c.id],cov[c.id,],100000)



bayeHR <- function(res,thin,id=10){
	qu <- round(quantile(res[thin,id],p=c(0.5,0.025,0.975)),2)
	equ <- round(exp(qu),2)
	
	r1 <- paste(qu[1]," (",qu[2]," - ",qu[3],")",sep="")
	r2 <- paste(equ[1]," (",equ[2]," - ",equ[3],")",sep="")
	c(r1,r2)
}

bayeHR(res.ln.neg)
bayeHR(res.ln.pos)
bayeHR(res.t.2)
bayeHR(res.t.3)
bayeHR(res.t.4)
bayeHR(res.g.1)
bayeHR(res.g.2)
bayeHR(res.g.3)
bayeHR(res.c.1)
bayeHR(res.c.2)


### comparison
rand_trial_cm <- coxph(s.ob~(factor(nodes)+factor(t))+factor(grade)+(lca199)+rand_arm_chk,data=espac4)

sg_ln.neg <- coxph(s.ob~(factor(t))+factor(grade)+(lca199)+rand_arm_chk,data=espac4[which(espac4$nodes==1),])

sg_ln.pos <- coxph(s.ob~(factor(t))+factor(grade)+(lca199)+rand_arm_chk,data=espac4[which(espac4$nodes==2),])

sg_t.2 <- coxph(s.ob~factor(nodes)+factor(grade)+(lca199)+rand_arm_chk,data=espac4[which(espac4$t==2),])

sg_t.3 <- coxph(s.ob~factor(nodes)+factor(grade)+(lca199)+rand_arm_chk,data=espac4[which(espac4$t==3),])

sg_t.4 <- coxph(s.ob~factor(nodes)+factor(grade)+(lca199)+rand_arm_chk,data=espac4[which(espac4$t==4),])

sg_g.1 <- coxph(s.ob~factor(nodes)+factor(t)+(lca199)+rand_arm_chk,data=espac4[which(espac4$grade==1),])

sg_g.2 <- coxph(s.ob~(factor(nodes)+factor(t))+(lca199)+rand_arm_chk,data=espac4[which(espac4$grade==2),])

sg_g.3 <- coxph(s.ob~(factor(nodes)+factor(t))+(lca199)+rand_arm_chk,data=espac4[which(espac4$grade==3),])

sg_c.1 <- coxph(s.ob~(factor(nodes)+factor(t)+factor(grade))+(lca199)+rand_arm_chk,data=espac4[which(espac4$lca199 < log(20)),])

sg_c.2 <- coxph(s.ob~(factor(nodes)+factor(t)+factor(grade))+(lca199)+rand_arm_chk,data=espac4[which(espac4$lca199 >= log(20)),])





summary(sg_ln.neg)
summary(sg_ln.pos)

summary(sg_t.2)
summary(sg_t.3)
summary(sg_t.4)

bayeHR(res.t.2)
bayeHR(res.t.3)
bayeHR(res.t.4)

summary(sg_g.1)
summary(sg_g.2)
summary(sg_g.3)

bayeHR(res.g.1)
bayeHR(res.g.2)
bayeHR(res.g.3)

summary(sg_c.1)
summary(sg_c.2)

bayeHR(res.c.1)
bayeHR(res.c.2)



rand_trial_cm

e4[1:3,]

##############################################
##############################################


####### Figure for Manuscript

### Cox model and median survival
setEPS()
postscript("modelfit.eps",height=10,width=10)
par(mfrow=c(2,2))
hist(espac3$lp,col=4,breaks=20,xlim=c(0,4),xlab="",main="Prognostic Indicators: Gem (ESPAC-3)")
abline(v=quantile(espac3$lp,p=c(0.16,0.5,0.84)),lwd=3,col="hotpink",lty=3)
text(0.25,60,c("a)"),cex=2,font=2)

hist(e4gcap$lp,col=2,breaks=20,xlab="",xlim=c(0,4),main="Prognostic Indicators: GemCap (ESPAC-4)")
abline(v=quantile(e4gem$lp,p=c(0.16,0.5,0.84)),lwd=3,col=5,lty=3)
text(0.25,27,c("b)"),cex=2,font=2)


plot(res[thin,10],lwd=1,typ="s",col=4,main="Trace of Posterior Distribution",ylim=c(-1.5,1))
text(500,.75,c("c)"),cex=2,font=2)

plot(density(res[thin,10],bw=0.025),lwd=3,col=4,main="Density of Posterior Distribution")
abline(v=-0.209,col=2,lwd=3)
abline(v=c(-0.55, 0.074),col=2,lwd=3,lty=2)
text(-1,2.5,c("d)"),cex=2,font=2)
dev.off()
#dev.copy2pdf(file="model_fit.pdf")



### Summarising results

qu <- round(quantile(exp(res[thin,10]),p=c(0.025,0.5,0.975)),2)

r1 <- paste(round(mean(res[thin,10]),3)," (",round(sd(res[thin,10]),3),")",sep="")
r2 <- paste(qu[2]," (",qu[1],"-",qu[3],")",sep="")

co <- summary(rand_trial_cm)$coefficients

###
exp(quantile(res[,11],p=c(0.025,0.5,0.975)))
summary(rand_trial_cm)

thin <- seq(5000,25000,by=5)

nrow(e4gcap)
layout(matrix(c(1,2,3,3),2,2))
plot(res[thin,11],typ="s",main="Trace of (log) HR",ylab="(log) HR",col="darkorchid")
acf(res[thin,11])


mn.lp <- mean(e4gcap$lp)
S <- S0^exp(mn.lp)

sfit <- survfit(s.ob~1,data= e4gcap);sfit
plot(sfit,conf.int=F,col=4,lwd=4,xlim=c(0,60),xlab="Time (Months)",ylab="Overall Survival")
lines(time,S,col=2,lwd=4,lty=2,typ="l",xlim=c(0,60))	
lines(time,S^exp(mnres),col=5,lwd=4,lty=2,typ="l",xlim=c(0,60))	
legend(20,0.95,c("gemPSC survival estimates","ESPAC-4 - GemCap KM estimates"),bty="n",col=c(4,2),lwd=4)


dev.copy2pdf(file="HRfig.pdf")

exp(median(res[,11]))
round(quantile(exp(res[thin,11]),p=c(0.5,0.025,0.975)),2)


#############################################










#### Dump

## comparison of linear predictors
e4gem$lp_3 <- mm.e4gem%*%coef(fpm)[-c(1:3)]


est1 <- round(coef(fpm),2)
se1 <- round(sqrt(diag(vcov(fpm))),2)
res1 <- paste(est1," (",se1,")",sep="")

est2 <- round(coef(e4gem.fpm),2)
se2 <- round(sqrt(diag(vcov(e4gem.fpm))),2)
res2 <- paste(est2," (",se2,")",sep="")

fpm.comp.res <- cbind(names(coef(fpm)),res1,res2)

write.csv(fpm.comp.res,"fpm.comp.res.csv")







