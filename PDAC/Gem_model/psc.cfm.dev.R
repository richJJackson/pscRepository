### Model Functions

setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()


### Getting Model and data
setwd("~/Documents/GitHub/pscLibrary/PDAC/Gem_Vs_GemCap")
load("flexParaGem.R")
e4 <- read.csv("e4_data_cohort.csv")

pscfit(fpm,e4)
CFM <- pscCFM(fpm)
pscfit(CFM,e4)


#### data comp

CFM$datavis

DC <- e4


vc <- visComp(CFM,DC)

ggarrange(plotlist=vc)

    
    
    
facVisComp <- function(p,x){
  
  old.d <- p$data;old.d
  tit <- p$labels$title
  old.d$source="CFM"
  dbnew <- data.frame(table(x));dbnew
  dbnew$source <- "DC"
  names(dbnew) <- names(old.d)
  df <- rbind(old.d,dbnew)
  
  cls <- brewer.pal(max(3,nrow(df)),"BuGn")
  
  
  p <- ggplot(data=df,aes(fill=x,values=Freq))+
    geom_waffle(color="white",size=0.33,n_rows=10)+
    theme_void()+
    scale_fill_manual(values = cls) +
    facet_wrap(~source,ncol=1)+
    ggtitle(tit)+
    theme(legend.position="right") +
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.07))
  p
}





    if(cls%in%c("numeric","integer")){
      
      
    }
    
    
    }
  
  
  lapply(cfmVis,attributes)  
  
  CFM$cov_class
}









bin.mod <- psc::bin.mod
DC <- psc::data

pscfit(bin.mod,DC)
CFM <- psc.cfm(bin.mod)
pscfit(CFM,DC)



#### 
CFM <- psc.cfm(bin.mod)
CFM$datavis$t


################
################
################
################


####### Data summary functions

install.packages("gtsummary")
library(gtsummary)
library(RColorBrewer)
library(waffle)

ramodel.frame(bin.mod)
model.frame(fpm)
cfm <- bin.mod

cfmDataSumm(bin.mod)
dv <- cfmDataVis(bin.mod)









############################################################################
############################################################################
############################################################################
############################################################################
############################################################################




######### data_match
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
  
  for(i in 1:length(cls)){
    
    con <- which(names(dc.data)%in%nm[i]);con
    if(length(con)==0) stop("DC missing covariate included in CFM")
    
    old <-dc.data[,con];old
    new <- old
    cl <- cls[[i]];cl
    
    if(cl%in%"character"){
      new <- factor(old);new
      warning(nm[i]," specified as a character in the model, consider respecifying 
                as a factor to ensure categories match between CFM and DC")
    }  
    
    if(cl%in%c("factor")){
      new <- factor(old);new
      if(!any(levels(new)%in%lev[[i]])) stop(paste("Factor levels in",nm[i],"not
                                                 represented in model"))
      att <- list("levels"=lev[[i]],class=cl)
      attributes(new) <- att
    }
    
    if(cl%in%c("numeric","integer")){
      new <- as.numeric(as.character(old))
    }
    
    dc.new[,con] <- new;dc.new
    rm(cl)
  }
  
  ret <- dc.new[,which(names(dc.new)%in%nm)]
  ret
}



########## ModelExtract
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


modelExtract.glm <- function(CFM){
  
  ### Model class
  mod_cls <- class(CFM)
  
  ## Model Ceofficients
  co <- CFM$coefficients;co
  fam <- CFM$family
  sig <- vcov(CFM)
  
  ### Data names
  form <- formula(CFM)
  mf <- model.frame(CFM)
  cls <-lapply(mf,class);cls
  
  char.id <- which(cls=="character")
  lev <- lapply(mf,levels);lev
  nm <- names(mf);nm
  
  ret <- list("mod_class"=mod_cls,"terms"=nm,"cov_class"=cls,"cov_lev"=lev,
              "cov_co"=co,"sig"=sig,"formula"=form,"family"=fam,"out.nm"=nm[1])
  ret
  
}





dataComb.glm <- function(CFM,DC,id=NULL,trt=NULL,cfmOb=F){
  
  ### removing response and weights
  if(!cfmOb) model_extract <- modelExtract(CFM);model_extract
  if(cfmOb) model_extract <- CFM
  
  ### Getting term names (and removing outcome and 'weights')
  term.nm <- model_extract$terms;term.nm
  out.nm <- term.nm[1]
  term.nm <- term.nm[-1];term.nm
  
  ### ERROR CHECK: Selecting data from DC
  data_unavail_id <- which(!term.nm %in% names(DC))
  data_unavail <- term.nm[data_unavail_id]
  if (length(data_unavail_id) != 0)
    stop(paste("Covariate '", data_unavail, "' is included in the model but not the dataset",
               sep = ""))
  out.id <- which(names(DC) %in% c(out.nm))
  if (length(out.id) != 1)
    stop(paste("Please ensure covariates for the outcome labelled",
               out.nm, "is included"))
  
  
  ### Adding treatment variable (if not null)
  if(!is.null(trt)) {
    if("trt"%in%names(DC)){
      DC <- DC[,-which(names(DC)=="trt")]
    }
    term.nm <- c(term.nm,"trt")
  }
  
  #### Selecting subgroup (if 'id' is specified)
  if(!is.null(id)){
    DC <- DC[id,]
    trt <- trt[id]
  }
  
  ### Removing missing data
  miss.cov <- which(is.na(DC),arr.ind=T)[,1]
  miss.trt <- which(is.na(trt))
  miss.id <- union(miss.cov,miss.trt)
  
  if(length(miss.id)>0) {
    DC <- DC[-miss.id,]
    trt <- trt[-miss.id]
    warning(paste(length(miss.id),"rows removed due to missing data in dataset"))
  }
  
  
  ## Matching data between DC and CFM
  cls <- model_extract$cov_class;cls
  lev <- model_extract$cov_lev;lev
  DCcov <- data_match(cls,lev,DC);DCcov[1:4,];trt[1:4]
  
  ## Defining outcome
  out <- data.frame(out.nm = DC[, which(names(DC) == out.nm)])
  names(out) <- out.nm
  dc_mm <- model.matrix(model_extract$formula,data=DCcov)
  
  ### Adding in 'trt' (if required)
  if(!is.null(trt)) dc_mm <- cbind(dc_mm,"trt"=DC$trt)
  
  ### returning results
  ret <- list("model.type"=class(CFM),"model_extract"=model_extract,"cov"=dc_mm,"outcome"=out)
  ret
  
}


pscEst.pscCFM <- function(CFM,DC_clean,nsim, start, start.se,trt=NULL){
  if("flexsurvreg"%in%CFM$mod_class) ret <- pscEst.flexsurvreg(CFM,DC_clean,nsim, start, start.se,trt=trt)
  if("glm"%in%CFM$mod_class) ret <- pscEst.glm(CFM,DC_clean,nsim, start, start.se,trt=trt)
  return(ret)
}



psc.cfm <- function(cfm,dataSumm=T,dataVis=T){
  me <- modelExtract(cfm)
  
  if(dataSumm){
    me$datasumm <- cfmDataSumm(cfm)
  }
  
  if(dataVis){
    me$datavis <- cfmDataVis(cfm)
  }
  
  class(me) <- c("pscCFM")
  return(me)
  me
}

