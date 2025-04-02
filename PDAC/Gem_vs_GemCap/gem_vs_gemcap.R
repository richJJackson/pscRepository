### Synthetic Controls - ESPAC3 and ESPAC4

#### library
library(psc)

setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()
library(knitr)
library(RColorBrewer)
library(ggplot2)
library(waffle)
library(ggpubr)


### Getting Model and data
setwd("~/Documents/GitHub/pscLibrary/PDAC/Gem_Vs_GemCap")
load("flexParaGem.R")

e4 <- read.csv("e4_data_cohort.csv")


#### Viewing package

#install.packages("pkgnet")
#library(pkgnet)
#report <- CreatePackageReport(pkg_name="psc",pkg_path="/Users/richardjackson/Documents/GitHub/psc")

###################################################################
###### Summarizing FPM model



attributes(fpm)

model_extract <- modelExtract(CFM);model_extract
mf <- model_extract$model.frame

mf <- mf[,-1]
mf <- mf[1:15,-5]


class(mf)
lapply(mf,class)
mf[,1] <- as.numeric(as.character(mf[,1]))
mf[,2] <- as.numeric(as.character(mf[,2]))
mf[,3] <- as.numeric(as.character(mf[,3]))
mf[,4] <- as.numeric(as.character(mf[,4]))



psc_espac <- pscfit(fpm,e4)
psc_espac
plot(psc_espac)

##### View Data from model




###################################################################

#### Sub group analysis

id<- which(e4$nodes==1)
psc_espac_n0 <- pscfit(fpm,e4,id=id)

id<- which(e4$nodes==2)
psc_espac_n1 <- pscfit(fpm,e4,id=id)


#### Multiple Treatment comparison

trt <- rbinom(nrow(e4),1,0.5);trt<-trt+1
psc_espac_mtc <- pscfit(fpm,e4,trt=trt)
psc_espac_mtc


### Warning if trt contains 0!?



plot(psc_espac_n0)
plot(psc_espac_n1)

CFM <- fpm
DC <- e4
nsim <- 5000
id <- NULL
trt <- trt

DC_clean



pscEst.flexsurvreg

pscfit <- function (CFM, DC, nsim = 5000, id = NULL, trt = NULL) {
  
  ### Cleaning Data
  DC_clean <- dataComb(CFM, DC, id=id, trt = trt)
  
  ### Starting Parameters
  init <- initParm(CFM = CFM, DC_clean = DC_clean, trt = trt)
  start<- init$par
  start.se <- sqrt(solve(init$hess))
  
  ### MCMC estimation### MhessianCMC estimation
  mcmc <- pscEst.flexsurvreg(CFM = CFM, DC_clean = DC_clean, nsim = nsim,
                 start = init$par, start.se=start.se,trt = trt)
  
  
  ####################################
  ####################################
  
  #pscEst.flexsurvreg<-   function(CFM,DC_clean,nsim,start,start.se,trt=trt){
    
    cov_co <- DC_clean$model_extract$cov_co;cov_co
    haz_co <- DC_clean$model_extract$haz_co;haz_co
    sig <- DC_clean$model_extract$sig
    lam <- DC_clean$model_extract$lam
    kn <- DC_clean$model_extract$kn
    time <- DC_clean$out$time;time
    cen <- DC_clean$out$cen
    cov <- DC_clean$cov
    est <- c(haz_co,cov_co);est
    
    trt.con <- is.null(trt)

    ####### Bayesian Estimation
    beta <- start
    parm <- matrix(NA,nsim,length(est)+length(beta)+1)
    parm[1,]<- c(est,beta,NA);parm[1,]
    
    s1 <- rmvnorm(1000,start,start.se)
    plot(s1[,1],s1[,2])
    
    ### multiplier for target distribution
    mult <- (5+5*start.se);mult
    
    ## Progress Bar
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    
    for(n in 2:nsim){
      
      ## progress bar
      setTxtProgressBar(pb, n)
      
      ### Drawing Samples
      cand <- rmvnorm(1,est,sig)
      cand.beta <- rnorm(length(beta),start,start.se*mult) #Check this bit
      parm[n,] <- c(cand,cand.beta,NA)
      
      ### partitioning covariates into baseline hazard and coefficients
      DC_cand <- DC_clean
      DC_cand$model_extract$cov_co <- cand[-c(1:length(haz_co))]
      DC_cand$model_extract$haz_co <- cand[1:length(haz_co)]
      
      
      ### cadidate log Hazard Ratios
      beta.old <- parm[n-1,-c(1:length(cand),ncol(parm))];beta.old
      beta.new <- parm[n,-c(1:length(cand),ncol(parm))];beta.new
      
      ### Prior contribution
      pr.cand <- -dmvnorm(cand,est,sig,log=T)
      
      pr.old <- dmvnorm(beta.old,rep(0,length(beta)),diag(length(beta))*1000,log=T)
      pr.new <- dmvnorm(beta.new,rep(0,length(beta)),diag(length(beta))*1000,log=T)
      
      ### Likelihood evaluation
      if(trt.con){
        l.old <- lik.flexsurvreg(beta.old,DC_cand) + pr.old + pr.cand
        l.new <- lik.flexsurvreg(beta.new,DC_cand)  + pr.new + pr.cand
      }
      
      if(!trt.con){
        l.old <- lik.flexsurvreg.mtc(beta.old,DC_cand) + pr.old + pr.cand;l.old
        l.new <- lik.flexsurvreg.mtc(beta.new,DC_cand)  + pr.new + pr.cand
      }
      
      ### Accept/Reject
      parm[n,ncol(parm)] <- l.new
      if(!acc(l.old,l.new)) {
        parm[n,-c(1:length(cand),ncol(parm))] <- beta.old
        parm[n,ncol(parm)] <- l.old
      }
      
    }
    
    parm
  }
  
  
  
  ####################################
  ####################################
  
  
  ### Formatting results
  covnm <- "beta"
  if (!is.null(trt)) {
    df <- data.frame(DC_clean$cov)
    ft <- factor(df$trt)
    covnm <- paste("beta", levels(ft), sep = "_")
  }
  
  ### Returning results 
  mcmc <- data.frame(mcmc)
  names(mcmc) <- c(colnames(DC_clean$model_extract$sig), covnm,
                   "DIC")
  psc.ob <- list(model.type = class(CFM), DC_clean = DC_clean,
                 posterior = mcmc)
  class(psc.ob) <- "psc"
  return(psc.ob)
}








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





