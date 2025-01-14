### Synthetic Controls - ESPAC3 and ESPAC4

#### library
library(psc)
library(knitr)


### Getting Model
setwd("~/Documents/GitHub/pscLibrary/PDAC/Model")
load("flexParaGem.R")
load("espac3fpm.R")


### Getting Data
setwd("~/Documents/GitHub/pscLibrary/PDAC/Data")
espac4 <- read.csv("espac4.csv")


###################################################################
espac3$t <- factor(espac3$t)
espac3$grade <- factor(espac3$grade)
espac3$nodes <- factor(espac3$nodes)

library(flexsurv)
fpm <- flexsurvspline(s.ob~t+grade+nodes+lca199,data=espac3,k=3)
save(fpm,file="flexParaGem.R")




###################################################################
###### Summarizing FPM model

fpm.res <- fpm$res
fpm.eres <- exp(fpm.res)

r1 <- paste(round(fpm.res[,1],2)," (",round(fpm.res[,4],2),")",sep="")
r2 <- paste(round(fpm.eres[,1],2)," (",round(fpm.eres[,2],2)," - ",round(fpm.eres[,3],2),")",sep="")

fpm.kab <- cbind(r1,r2)
colnames(fpm.kab) <- c("Est (se)","HR (95% CI)")
rownames(fpm.kab) <- rownames(fpm.res)
kable(fpm.kab)

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




e4gem[1:3,]

psc_espac <- pscfit(fpm,e4gcap)
psc_espac
exp(-0.239)
p <- plot(psc_espac)

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





