### Model Functions

setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()


### Getting Model and data
setwd("~/Documents/GitHub/pscLibrary/PDAC/Gem_Vs_GemCap")
load("flexParaGem.R")

e4 <- read.csv("e4_data_cohort.csv")

e4[1:3,]

### psc.cfm OBJECT

### INCLUDE - Model Extract


me <- modelExtract(fpm)
dataComb.flexsurvreg(fpm,e4)


pscfit <- function (CFM, DC, nsim = 5000, id = NULL, trt = NULL) {
  
  ### Cleaning Data
  DC_clean <- dataComb(CFM, DC, id=id, trt = trt)
  
  ### Starting Parameters
  init <- initParm(CFM = CFM, DC_clean = DC_clean, trt = trt)
  start<- init$par
  start.se <- sqrt(solve(init$hess))
  
  ### MCMC estimation### MhessianCMC estimation
  mcmc <- pscEst(CFM = CFM, DC_clean = DC_clean, nsim = nsim,
                 start = init$par, start.se=start.se,trt = trt)
  
  ### Formatting results
  covnm <- "beta"
  if (!is.null(trt)) {
    df <- data.frame(DC_clean$cov)
    ft <- factor(df$trt)
    covnm <- paste("beta", levels(ft), sep = "_")
  }
  
  mcmc <- data.frame(mcmc)
  names(mcmc) <- c(colnames(DC_clean$model_extract$sig), covnm,
                   "DIC")
  psc.ob <- list(model.type = class(CFM), DC_clean = DC_clean,
                 posterior = mcmc)
  class(psc.ob) <- "psc"
  return(psc.ob)
}


psc.cfm <- function(mod){
  
  me <- modelExtract(mod)
  class(me) <- c("pscCFM")
  return(me)
}

CFM <- psc.cfm(fpm)




DC <- e4


dataComb(CFM,DC,id=NULL,trt=NULL)





dataComb.pscCFM <- function(CFM,DC,id=NULL,trt=NULL,cfmOb=T){
  
  if("flexsurvreg"%in%CFM$mod_class) ret <- dataComb.flexsurvreg(CFM=CFM,DC=DC,id=NULL,trt=NULL,cfmOb=T)
  if("glm"%in%CFM$mod_class) ret <- dataComb.glm(CFM=CFM,DC=DC,id=NULL,trt=NULL,cfmOb=T)
  
  return(ret)
  }



modelExtract.flexsurvreg <- function(CFM){
  
  ### Model class
  mod_cls <- class(CFM)
  
  ## Model Ceofficients
  co <- CFM$coefficients
  k <- CFM$k
  kn <- CFM$knots
  sig <- vcov(CFM)
  lam <- (max(kn)-kn)/(max(kn)-min(kn))
  
  ### Data names
  form <- formula(CFM)
  mf <- CFM$data$m
  lev <- lapply(mf,levels)
  cls <-lapply(mf,class)
  nm <- names(mf)

  
  ### Cleaning model parameters
  n_haz_co <- k+2
  haz_co <- co[1:n_haz_co]
  cov_co <- co[(n_haz_co+1):length(co)]
  
  ret <- list("mod_class"=mod_cls,"terms"=nm,"cov_class"=cls,"cov_lev"=lev,
              "cov_co"=cov_co,"sig"=sig,"haz_co"=haz_co,"k"=k,"kn"=kn,
              "lam"=lam,"formula"=form)
  return(ret)
}



dataComb.flexsurvreg <- function(CFM,DC,id=NULL,trt=NULL,cfmOb=F){

  ### removing response and weights
  if(!cfmOb) model_extract <- modelExtract(CFM);model_extract
  if(cfmOb) model_extract <- CFM
  
  ### Getting term names (and removing outcome and 'weights')
  term.nm <- model_extract$terms
  term.nm <- term.nm[-c(1,length(term.nm))];term.nm

  ### ERROR CHECK: Selecting data from DC
  data_unavail_id  <- which(!term.nm%in%names(DC))
  data_unavail <- term.nm[data_unavail_id]
  if(length(data_unavail_id)!=0) stop(paste("Covariate '",data_unavail,"' is included in the model but not the dataset",sep=""))

  ### Making sure 'time' and 'cen' are present
  out.nm.trap <- which(names(DC)%in%c("time","cen"))
  if(length(out.nm.trap)!=2) stop("outcome covariates in data cohort should be named 'time' and 'cen'")
  
  ## Creating model matrix (adding resposne to covariates)
  DC <- data.frame(cbind(0,DC))
  names(DC)[1] <- names(mf)[1]
  DC[,1] <- Surv(DC$time,DC$cen)

  ### Adding treatment variable (if not null)
  if(!is.null(trt)) {
    if("trt"%in%names(DC)){
      DC <- DC[,-which(names(DC)=="trt")]
    }
    term.nm <- c(term.nm,"trt")
  }
  
  ## Matching data between DC and CFM
  cls <- model_extract$cov_class
  lev <- model_extract$cov_lev
  DCcov <- data_match(cls,lev,DC);DCcov[1:4,];trt[1:4]

  #### Selecting subgroup (if 'id' is specified)
  if(!is.null(id)){
    DCcov <- DCcov[id,]
    trt <- trt[id]
  }

  ### Removing missing data
  miss.cov <- which(is.na(DCcov),arr.ind=T)[,1]
  miss.trt <- which(is.na(trt))
  miss.id <- union(miss.cov,miss.trt)

  if(length(miss.id)>0) {
    DCcov <- DCcov[-miss.id,]
    trt <- trt[-miss.id]
    warning(paste(length(miss.id),"rows removed due to missing data in dataset"))
  }

  ### Creating model matrix based on new dataset
  out <- DCcov[,1]
  dc_mm <- model.matrix(model_extract$formula,data=DCcov)[,-1]
  
  ### Adding in 'trt' (if required)
  if(!is.null(trt)) dc_mm <- cbind(dc_mm,"trt"=trt)
  ret <- list("model_extract"=model_extract,"cov"=dc_mm,"outcome"=out)
  ret
  
}




data_match <- function(cls,lev,dc.data){
  
  ### duplicated namse in dc.data
  dup <- duplicated(names(dc.data))
  if(any(dup)) dc.data <- dc.data[,-which(dup)]
  
  ## Getting term names
  nm <- names(cls);nm
  
  ## removing 'weights' column if there
  wid <- which(nm=="(weights)");wid
  if(length(wid)>0) cls <- cls[-wid]
  
  # creating output
  dc.new <- dc.data
  
  for(i in 1: length(cls)){
    
    con <- which(names(dc.data)%in%nm[i]);con
    if(length(con)==0) stop("DC missing covariate included in CFM")
    
    old <-dc.data[,con];old
    new <- old
    cl <- cls[[i]];cl
    
    if(cl%in%c("character","factor")){
      new <- factor(old);new
      if(!any(levels(new)%in%lev[[i]])) stop(paste("Factor levels in",nm[i],"not
                                                 represented in model"))
      att <- list("levels"=lev[[i]],class="cl")
      attributes(new) <- att
    }
    
    if(cl%in%c("numeric","integer")){
      new <- as.numeric(as.character(old))
    }
    
    dc.new[,con] <- new
    rm(cl)
  }
  
  ret <- dc.new[,which(names(dc.new)%in%nm)]
  ret
}

