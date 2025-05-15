


###### Validation of Gemcitabine Model


########### Function used to validate fpm model

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


#########

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



