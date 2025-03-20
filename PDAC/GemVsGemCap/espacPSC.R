### Synthetic Controls - ESPAC3 and ESPAC4

#### library
library(psc)
setwd("/Users/richardjackson/Documents/GitHub/pscVis/pscVis")
devtools::load_all()
setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()
library(knitr)
library(RColorBrewer)
library(ggplot2)
library(waffle)
library(ggpubr)
library(table1)
setwd("/Users/richardjackson/Documents/GitHub/toolsRJJ")
devtools::load_all()

### Getting Model and data
setwd("~/Documents/GitHub/pscLibrary/PDAC/GemVsGemCap")
load("flexParaGem.R")

e4 <- read.csv("e4_data_cohort.csv")



out <- data.frame(fpm$data$Y)
survfit(Surv(out$time,out$status)~1)

###################################################################
###### Summarizing FPM model

fpm.res <- fpm$res
fpm.eres <- exp(fpm.res)


### Summarising data
names(fpm$data$m)[2:5] <- c("T Stage","Tumor Grade","Lymph Node","(log) CA19.9")
data <- fpm$data$m[,2:5]
data$ap <- "All Pts."

stab <- summaryTable(data[,1:4],by=data$ap)[,-4]
rownames(stab) <- NULL

kable(stab,format="html")


#### Summary of DaNULL#### Summary of Dataset
png("e3_data.png", width=600, height=1200, bg = "transparent"); 
ggarrange(vis[[1]],vis[[2]],vis[[3]],vis[[4]],nrow=4)
dev.off()

### Summary Table
r1 <- paste(round(fpm.res[,1],2)," (",round(fpm.res[,4],2),")",sep="")
r2 <- paste(round(fpm.eres[,1],2)," (",round(fpm.eres[,2],2)," - ",round(fpm.eres[,3],2),")",sep="")

fpm.kab <- cbind(r1,r2)
colnames(fpm.kab) <- c("Est (se)","HR (95% CI)")
rownames(fpm.kab) <- rownames(fpm.res)
kable(fpm.kab,format="html")

### Summary Figure
png("e3_gem_fpm.png", width=600, height=600, bg = "transparent"); 
par(bty="l",mar=c(5,5,2,1))
plot(fpm,col=5,lwd=4,xlim=c(0,60),lwd.obs=3,col.obs="darkorchid",cex.axis=1.4,
     cex.lab=1.5,ylab="Overall Survival",xlab="Time (Months)",font.lab=3)
abline(h=seq(0,1,by=0.1),v=seq(0,60,by=12),lty=2,col="gray")
legend(40,0.9,c("Data","Model"),col=c("darkorchid",5),lwd=4,bty="n",cex=1.5)
dev.off()




###### Validation

val <- cfmValid.flexsurvreg(fpm)


png("e3_gem_discr.png", width=600, height=600, bg = "transparent"); 
plot(val$discrim$sfit,col=c(4,5,6,"darkorchid"),lwd=4,xlim=c(0,60),lwd.obs=3,col.obs="darkorchid",cex.axis=1.4,
     cex.lab=1.5,ylab="Overall Survival",xlab="Time (Months)",font.lab=3)
abline(h=seq(0,1,by=0.1),v=seq(0,60,by=12),lty=2,col="gray")
legend(40,0.9,c("Risk Group 1","Risk Group 2","Risk Group 3","Risk Group 4"),col=c(4,5,6,"darkorchid"),lwd=4,bty="n",cex=1.5)
dev.off()


co <- summary(val$discrim$cm)$coefficients
lo <- co[,1]-1.96*co[,3]
up <- co[,1]+1.96*co[,3]

est <- round(co[,1],3)
hr <- round(exp(est),3)
se <- round(co[,3],3)
lo <- round(exp(lo),3)
up <- round(exp(up),3)


r1 <- paste(est," (",se,")",sep="")
r2 <- paste(hr," (",lo,", ",up,")",sep="")

tb <- rbind("",cbind(r1,r2))
rownames(tb) <- c("Risk Group 1","Risk Group 2","Risk Group 3","Risk Group 4")
colnames(tb) <- c("est (se)","HR (95% CI)")

kable(tb,format="html")



calib_c <- round(val$calib$c,3)
calib_slope <- round(summary(val$calib$slope)$coefficients,3)


calib_c <- paste(calib_c[1]," (",calib_c[2],")",sep="")
calib_slope <- paste(calib_slope[1]," (",calib_slope[3],")",sep="")

tb2 <- rbind(calib_c,calib_slope)
rownames(tb2) <- c("C-Statistic","Slope")
colnames(tb2) <- "est (se)"
kable(tb2,format="html")



attributes(val)

cfmValid.flexsurvreg <- function(cfm,exData=NULL){
  
  ret <- list()
  ret$detail <- list
  ret$discrim <- list()
  ret$calib <- list()
  
  ## Detail
  ret$detail <- "Internal validation based on data used to train model"
  
  if(!is.null(exData)){
    ret$detial <- "External validation based on dataset other than that used to train the model"
  }
  
  
  ## Discrimination
  me <- modelExtract.flexsurvreg(cfm)
  cov <- model.matrix(cfm)
  
  ###
  ##
  # Add in external data here
  ##
  ###
  
  lp <- t(me$cov_co%*%t(cov))
  
  ## Creating groups
  pred_quant <- quantile(lp,c(0.15,0.5,0.85))
  pred_grp <- cut(lp,c(-Inf,pred_quant,Inf),labels=c("risk_grp_1","risk_grp_2","risk_grp_3","risk_grp_4"))
  
  ##
  out <- cfm$data$m[,1]
  sfit <- survfit(out~pred_grp)
  cm <- coxph(out~pred_grp)
  
  ret$discrim$sfit <- sfit
  ret$discrim$cm <- cm
  
  ## Calibration
  
  # concordance & slope
  c <- summary(cm)$concordance
  lp_slope <- coxph(out~lp)
  
  ret$calib$c <- c
  ret$calib$slope <- lp_slope
  
  ## returning object
  ret
}




###################################################################


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
espac4$lca199 <- log(espac4$ca199_post+1)

e4 <- espac4[,which(names(espac4)%in%c("stime","event","nodes","grade","t","lca199"))]
names(e4) <- c("cen","grade","time","nodes","lca199","t")

e4gem <- e4[which(espac4$rand_arm_chk==1),]
e4gcap <- e4[which(espac4$rand_arm_chk==2),]

### removing missing data
e4gem <- e4gem[-unique(which(is.na(e4gem),arr.ind=T)[,1]),]
e4gcap <- e4gcap[-unique(which(is.na(e4gcap),arr.ind=T)[,1]),]



psc_espac <- pscfit(fpm,e4gcap)
psc_espac
exp(-0.239)
p <- plot(psc_espac)

###################################################################



id<- which(e4gcap$nodes==1)
psc_espac_n0 <- pscfit(fpm,e4gcap,id=id)
sd(psc_espac_n0$posterior$beta)
exp(coef(psc_espac_n0))

id<- which(e4gcap$nodes==2)
psc_espac_n1 <- pscfit(fpm,e4gcap,id=id)
psc_espac_n1
sd(psc_espac_n1$posterior$beta)
exp(coef(psc_espac_n1))


table(e4gcap$grade)

id<- which(e4gcap$t==2)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))

id<- which(e4gcap$t==3)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))

id<- which(e4gcap$t==4)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))


id<- which(e4gcap$grade==2)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))

id<- which(e4gcap$grade==3)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))


e4gcap[1:3,]

id<- which(exp(e4gcap$lca199)<20)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))

id<- which(exp(e4gcap$lca199)>=20)
psc_tmp<- pscfit(fpm,e4gcap,id=id)
psc_tmp
sd(psc_tmp$posterior$beta)
exp(coef(psc_tmp))
###################################################################

espac4[1:3,]




summary(coxph(s.ob~r_status+rand_arm_chk,data=espac4[id,]))


id <- which(espac4$ca199_post<20)
  
plot(psc_espac$posterior$beta)
sd(psc_espac$posterior$beta)

library(survminer)
p <- plot.psc.flexsurvreg(psc_espac)
ggsave(p,file="gg_plot.png")


plot.psc.flexsurvreg <- function(x, ...){
  
  # Binding local varaibles
  S <- trt <- NULL
  
  med <- coef(x)
  med <- med[-nrow(med),1]
  
  ## defining treatment (for multiple treatment comparisons)
  mtc.cond <- "trt"%in%colnames(x$DC_clean$cov)
  trt <- rep(1,nrow(x$DC_clean$cov))
  if(mtc.cond) trt <- factor(x$DC_clean$cov[,which(colnames(x$DC_clean$cov)=="trt")])
  
  
  ### Getting model survival estimate
  s_fpm <- surv_fpm(x$DC_clean)
  s_data <- data.frame("time"=s_fpm$time,"S"=s_fpm$S)
  
  # plot
  out <- x$DC_clean$outcome
  out$trt <- trt
  sfit <- survfit(Surv(time,cen)~trt,data=out)
  sfit_plot <- ggsurvplot(sfit,data=out,legend="none",color="blue")$plot
  sfit_plot + 
    geom_line(data=s_data, aes(time,S,color="purple"),lwd=1.5) +
    labs(y="Survival Probability",x="Time2",color="Legend") +
    scale_colour_manual(name="Error Bars",values=cols)
  
}

plot.psc.flexsurvreg(psc_espac)
library(ggplot2)


###  Synthetic Control - carful not to set the sims too high - it takes a while!
res <- synthComp(fpm,tim,cen,cov,100000)

thin <- seq(5000,100000,by=10)
quantile(res[thin,10],p=c(0.5,0.025,0.975))
exp(quantile(res[thin,10],p=c(0.5,0.025,0.975)))





