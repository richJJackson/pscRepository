###################################################################
###################################################################
### Synthetic Controls - ESPAC3 and ESPAC4

#### library
library(psc)


rm(list=ls())
setwd("/Users/richardjackson/Documents/GitHub/psc")
devtools::load_all()
library(knitr)
library(RColorBrewer)
library(ggplot2)
library(waffle)
library(ggpubr)
library(gtsummary)
library(forestplot)

###################################################################
###################################################################


### Getting Model
setwd("~/Documents/GitHub/pscLibrary/PDAC/Gem_model")
load("CFM.Rds")


## View summary of data
ds <- CFM$datasumm$summ_Table;ds
gt::gtsave(as_gt(ds),"dataSumm.pdf")

## plot data 
png("gem_model_dv.png")
ggarrange(plotlist=CFM$datavis,nrow=2,ncol=2)
dev.off()


###################################################################
###################################################################


### Getting Data
setwd("~/Documents/GitHub/pscLibrary/PDAC/Gem_Vs_GemCap")
e4 <- read.csv("e4_data_cohort.csv")


### Summary of Data
e4.cov <- e4[,-c(1:2)]
e4_summ <- tbl_summary(e4.cov,missing="ifany")
gt::gtsave(as_gt(e4_summ),"e4DataSumm.pdf")


### Visual Comparison DC (e4) against the CFM
vc <- visComp(CFM,e4)
ggarrange(plotlist=vc,ncol=2,nrow=2)


model_extract <- CFM

############# Overall Comparison
psc_all <- pscfit(CFM,e4)

### Plotting overall effect
plot(psc_all)

coef(psc_all)



############## Sub-Group analysis


### Nodes
id_n1 <- which(e4$nodes==1)
psc_n1 <- pscfit(CFM,e4,id=id_n1)

id_n2 <- which(e4$nodes==2);id_n2
psc_n2 <- pscfit(CFM,e4,id=id_n2)

pn1 <- plot(psc_n1)
pn2 <- plot(psc_n2)

ggarrange(pn1,pn2)



### grade
id_g1 <- which(e4$grade==1)
psc_g1 <- pscfit(CFM,e4,id=id_g1)

id_g2 <- which(e4$grade==2);id_g2
psc_g2 <- pscfit(CFM,e4,id=id_g2)

id_g3 <- which(e4$grade==3);id_g3
psc_g3 <- pscfit(CFM,e4,id=id_g3)

pg1 <- plot(psc_g1)
pg2 <- plot(psc_g2)
pg3 <- plot(psc_g3)

ggarrange(pg1,pg2,pg3)



### t
id_t2 <- which(e4$t==2)
psc_t2 <- pscfit(CFM,e4,id=id_t2)

id_t3 <- which(e4$t==3);id_t3
psc_t3 <- pscfit(CFM,e4,id=id_t3)

id_t4 <- which(e4$t==4);id_t4
psc_t4 <- pscfit(CFM,e4,id=id_t4)

pt2 <- plot(psc_t2)
pt3 <- plot(psc_t3)
pt4 <- plot(psc_t4)

ggarrange(pt2,pt3,pt4)


e4[1:3,]

### lca199
id_c1 <- which(e4$lca199<2.5)
psc_c1 <- pscfit(CFM,e4,id=id_c1)

id_c2 <- which(e4$lca199>=2.5)
psc_c2 <- pscfit(CFM,e4,id=id_c2)

pc1 <- plot(psc_c1)
pc2 <- plot(psc_c2)

ggarrange(pc1,pc2)

#############################
### combined Results
res <- data.frame(rbind(
  coef(psc_n1)[1,],
  coef(psc_n2)[1,],
  coef(psc_g1)[1,],
  coef(psc_g2)[1,],
  coef(psc_g3)[1,],
  coef(psc_t2)[1,],
  coef(psc_t3)[1,],
  coef(psc_t4)[1,],
  coef(psc_c1)[1,],
  coef(psc_c2)[1,],
  coef(psc_all)[1,]
))


names(res) <- c("median","lower","upper","p0","p1")
res$label <- c("N1","N2","grade1","grade2","grade3","t1","t2","t3","c1","c2","All Pts.")

names(res)

### Forrest Plot


res |>
  forestplot(mean=median,labeltext = label,is.summary=c(F,F,F,F,F,F,F,F,F,F,T)) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue")|>
  fp_set_zebra_style("#EFEFEF")













