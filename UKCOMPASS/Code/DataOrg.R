### Sorting data for analysis

#Saving 'Home'
hm <- getwd()

## Setting working directory for data
setwd("Data")

### laoding data
load("ukc_ri.r")

## Defining outcome
which(!is.na(ri$late_opDt))


ri_late$opTm <- as.numeric(as.Date(ri_late$late_opDt,"%Y-%m-%d")-as.Date(ri_late$T0,"%d%b%Y"))/30.44
ri_late$out <- ((!is.na(ri_late$opTm))|ri_late$cen==1)
#ri_late$out <- as.numeric(ri_late$opTm<ri_late$stime)
ri_late$out[which(is.na(ri_late$out))] <- 0
ri_late$t2 <- factor(ri_late$tx,labels=c("OR","EV","EV"))

ri_late$I_sex <- factor(ri_late$I_sex)
ri_late$asa_cat <- factor(ri_late$asa_cat)

data <- ri_late
save(data,file="analData.R")



