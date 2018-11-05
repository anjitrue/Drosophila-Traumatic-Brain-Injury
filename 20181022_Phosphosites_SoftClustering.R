#install.packages("dplyr")
#install.packages("e1071")
#install.packages("foreign")
#install.packages("MASS")
#install.packages("Matrix")
#install.packages("mgcv")
#install.packages("survival")

library(BiocManager)
BiocManager::install("Mfuzz", version = "3.8")

library(e1071)
library(dplyr)
library(Mfuzz)


# Data In -----------------------------------------------------------------

setwd("P:\\GMW_20171110_Hornberger_Human\\No Isoform Search\\Code")

data <- read.csv("20181021 Subject Normalized Merged Data.csv")

colnames(data)

exprs <- data[,31:35]

rownames(exprs) <- paste(data[,2],data[,4])

as.matrix(exprs)

# mFuzz Soft Clustering----------------------------------------------------

require(Biobase)
exprValues <- new("ExpressionSet", exprs = as.matrix(exprs))

exprValues.r <- filter.NA(exprValues, thres = 0.9)

exprValues.f = fill.NA(exprValues.r, mode="mean")

tmp = filter.std(exprValues.f, min.std = 0.1)

exprValues.s = standardise(tmp)

m1 <- mestimate(exprValues.s)

error <- NA
for(i in 2:15){
  c1 <- mfuzz(exprValues.s, c=i, m=m1)
  error <- rbind(error, c(i,c1$withinerror))
}
plot(error[,1], error[,2])

c1 <- mfuzz(exprValues.s, c=6, m=m1)

pdf(file = "20181022_Phosphosites_SoftClusters.pdf")
par(mar=c(3,3,1,1))
mfuzz.plot2(exprValues.s, cl=c1, mfrow = c(3,2), col.lab="black", x11=F)
dev.off()

a <- acore(exprValues.s, c1, min.acore=0.0)

cluster1 <- data.frame(a[[1]])
cluster2 <- data.frame(a[[2]])
cluster3 <- data.frame(a[[3]])
cluster4 <- data.frame(a[[4]])
#For Phospho Data
cluster5 <- data.frame(a[[5]])
cluster6 <- data.frame(a[[6]])

write.csv(cluster1, "cluster1Info.csv")
write.csv(cluster2, "cluster2Info.csv")
write.csv(cluster3, "cluster3Info.csv")
write.csv(cluster4, "cluster4Info.csv")
#For Phospho Data
write.csv(cluster5, "cluster5Info.csv")
write.csv(cluster6, "cluster6Info.csv")
