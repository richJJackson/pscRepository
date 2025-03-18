### output for webpage

### LOADING functions/libraries
library(readxl)
library(RColorBrewer)
library(waffle)
library(ggplot2)
library(ggpubr)

setwd("~/Documents/GitHub/psc")
devtools::load_all()
setwd("~/Documents/GitHub/pscVis/pscVis")
devtools::load_all()

### Getting TACE model
setwd("/Users/richardjackson/Documents/GitHub/pscLibrary/TACE_SIRT/Data")
load("1knotmodel.Rdata")
flexspline_model1



setwd("/Users/richardjackson/Documents/GitHub/pscLibrary/TACE_SIRT")

model.ob <- pscCFM(flexspline_model1)


ggarrange(model.ob$setting$vis[[1]],
          model.ob$setting$vis[[2]],
          model.ob$setting$vis[[3]],
          model.ob$setting$vis[[4]],nrow=4)
ggsave(file="TACE_dataPlot.png",width=6,height=12)          



summary(flexspline_model1)$coefficients
flexspline_model1$coefficients

png(filename="TACE_km1.png")
plot(flexspline_model1,ylab="Overall Survival (%)",xlab="Time (Months)",
     col=c(4,5),lwd=4,cex.lab=1.4,font.lab=3,cex.axis=1.3,cex.main=1.4,
     font.main=2,main="TACE Survival")
dev.off()


est <- round(c(model.ob$model$haz_co,model.ob$model$cov_co),3)
se <- round(sqrt(diag(model.ob$model$sig)),3)

tb <- cbind(est,se)

tb

library(knitr)
kable(tb,format="html")
ggarrange(model.ob)

tb1
