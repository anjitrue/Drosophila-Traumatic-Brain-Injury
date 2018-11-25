##################################Load Packages #################################################
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")

devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes=TRUE)

library(readr)
library(plyr)
library(dplyr)
library(pheatmap)
library(pcaMethods)
library(ggplot2)
library(devtools)
library(e1071)
library(dplyr)
library(Mfuzz)

##################################### Load Data ####################################################

#set working Directory
setwd("E:/Projects/Proteomics/DorsophilaHead_Experiment/")

#load data
proteinGroups_dros_hemo <- read_delim("E:/Projects/Proteomics/DorsophilaHead_Experiment/txt_hemo_plusFrac/proteinGroups.txt","\t", escape_double = FALSE, trim_ws = TRUE)
# dim(proteinGroups_dros_hemo) [1] 9746  131
#proteinGroups_dros_hemo <- data.frame(proteinGroups_dros_hemo)
proteinGroups_dros_heads <- read_delim("E:/Projects/Proteomics/DorsophilaHead_Experiment/txt_dros_heads_plusFrac/proteinGroups__dros_heads_plusFrac.txt","\t", escape_double = FALSE, trim_ws = TRUE)
#[1] 10290   338
##################################### Functions ####################################################

RemoveCotaminats <- function(x){
 x <- x[-grep("*",x$`Potential contaminant`),]
 x <- x[-grep("*",x$Reverse),]
 x <- x[-grep("*",x$`Only identified by site`),]
 x <- x[-grep("CON",x$`Protein IDs`),]
 x <- x[-grep("REV",x$`Protein IDs`),]
 x <- x[,-which(names(x) %in% c("Potential contaminant","Reverse","Only identified by site"))]
 return(x)
 
}

subsetLFQ <- function(x){
  y <- x[,which(names(x) %in% c("Protein IDs", "id","Gene names"))]
  z <- x[,grep("LFQ intensity",names(x))]
  z[z == 0] <- NA
  x <- data.frame(y,z)
  return(x)
}

remove.features.50percentcuttoff <- function (x) {
  features.missing = rowMeans(is.na(x)) 
  print(sum(features.missing > 0.50)) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  keep.features = which(features.missing <= 0.50) 
  print(paste0("Remaining Features: ", length(keep.features)))
  names(keep.features) = keep.features 
  
  remove.features = which(features.missing > 0.50)
  print(paste0("Features Removed: ", length(remove.features)))
  names(remove.features) = remove.features
  
  filtered = x[-which(rownames(x) %in% remove.features),]
  return(filtered)
}

remove.features.25percentcuttoff <- function (x) {
  features.missing = rowMeans(is.na(x)) 
  print(sum(features.missing > 0.75)) 
  features.missing.75more = rownames(x)[features.missing > 0.75] 
  
  keep.features = which(features.missing <= 0.25) 
  print(paste0("Remaining Features: ", length(keep.features)))
  names(keep.features) = keep.features 
  
  remove.features = which(features.missing > 0.75)
  print(paste0("Features Removed: ", length(remove.features)))
  names(remove.features) = remove.features
  
  filtered = x[-which(rownames(x) %in% remove.features),]
  return(filtered)
}

remove.features.10percentcuttoff <- function (x) {
  features.missing = rowMeans(is.na(x)) 
  print(sum(features.missing > 0.90)) 
  features.missing.10more = rownames(x)[features.missing > 0.90] 
  
  keep.features = which(features.missing <= 0.10) 
  print(paste0("Remaining Features: ", length(keep.features)))
  names(keep.features) = keep.features 
  
  remove.features = which(features.missing > 0.90)
  print(paste0("Features Removed: ", length(remove.features)))
  names(remove.features) = remove.features
  
  filtered = x[-which(rownames(x) %in% remove.features),]
  return(filtered)
}

fuzzyprep_imputation_included <- function(z)
{
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  # exclude proteins that have more than 50% of measurements missing
  exprValues.r <- filter.NA(exprValues, thres = 0.50)
  exprValues.f = fill.NA(exprValues.r, mode="median")
  tmp = filter.std(exprValues.f, min.std = 0.1)
  exprValues.s = standardise(tmp)
  return(exprValues.s)
}

fuzzyprep_usepreviousImputation <- function(z)
{
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  # exclude proteins that have more than 50% of measurements missing
  exprValues.r <- filter.NA(exprValues, thres = 0.50)
  tmp = filter.std(exprValues.r, min.std = 0.1)
  exprValues.s = standardise(tmp)
  return(exprValues.s)
}

##################################### Hemo Data ####################################################

proteinGroups_dros_hemo <- RemoveCotaminats(proteinGroups_dros_hemo) # dim(proteinGroups_dros_hemo) [1] 9181  128
proteinGroups_dros_hemo <- subsetLFQ(proteinGroups_dros_hemo) # dim(proteinGroups_dros_hemo) [1] 9181  12


sum(is.na(proteinGroups_dros_hemo)) #31443 missing values

filtered_dros_hemo_50percent <- remove.features.50percentcuttoff(proteinGroups_dros_hemo)
#[1] 3151
#[1] "Remaining Features: 6030"
#[1] "Features Removed: 3151"

filtered_dros_hemo_25percent <- remove.features.25percentcuttoff(proteinGroups_dros_hemo)
#[1] 3
#[1] "Remaining Features: 5571"
#[1] "Features Removed: 3"


filtered_dros_hemo_10percent <- remove.features.10percentcuttoff(proteinGroups_dros_hemo)
#[1] 0
#[1] "Remaining Features: 9181"
#[1] "Features Removed: 0"

even_index <- seq(2,8,2)
odd_index <- seq(1,8,2)

# set up for baysian pca imputation
hemo_50percent <- as.matrix(filtered_dros_hemo_50percent[,c(4:11)])
hemo_50percent_log2 <- log2(hemo_50percent)
rownames(hemo_50percent_log2) <- filtered_dros_hemo_50percent$id

pc_hemolog2 <- pca(hemo_50percent_log2, nPcs = 3, method = "bpca")
imputed_hemo <- completeObs(pc_hemolog2)

imputed_hemo_even <- imputed_hemo[,even_index]
imputed_hemo_odd <- imputed_hemo[,odd_index]


#subset for other modes of imputation
z_hemo <- proteinGroups_dros_hemo[,4:11]
z_hemo <- log2(z_hemo)
rownames(z_hemo) <- proteinGroups_dros_hemo$id
Z_hemo_even <- z_hemo[,even_index]
z_hemo_odd <- z_hemo[,odd_index]

# Hemo mFuzz Soft Clustering----------------------------------------------------
z = imputed_hemo_odd #change according to data input


exprValues.s_imput <- fuzzyprep_imputation_included(Z_hemo_even)
exprValues.s_bayesian <- fuzzyprep_usepreviousImputation(imputed_hemo_even)

m1 <- mestimate(exprValues.s_imput)
error <- NA

for(i in 2:15){
  c1 <- mfuzz(exprValues.s_imput, c=i, m=m1)
  error <- rbind(error, c(i,c1$withinerror))
}
plot(error[,1], error[,2])


#exprValues.s.nominstd = standardise(exprValues.r) #use data that is non std

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



##################################### Head Data ####################################################
proteinGroups_dros_heads <- RemoveCotaminats(proteinGroups_dros_heads) # [1] 9365  335
proteinGroups_dros_heads <- subsetLFQ(proteinGroups_dros_heads) # [1] 9365   38

z = as.matrix(proteinGroups_dros_heads[,5:length(proteinGroups_dros_heads)])

#order the columns
colnames(z) <- sub("LFQ.intensity.Heads_", "", colnames(z))
colnames(z) <- as.numeric(colnames(z))
z <- z[,order(as.integer(colnames(z)))]
z <- log2(z)

even_index <- seq(2,34,2)
odd_index <- seq(1,34,2)

z <- z[,odd_index]

sum(is.na(z)) #80327

filtered_dros_head_50percent <- remove.features.50percentcuttoff(proteinGroups_dros_heads)
#[1] 2190
#[1] "Remaining Features: 7129"
#[1] "Features Removed: 2190"

filtered_dros_head_25percent <- remove.features.25percentcuttoff(proteinGroups_dros_heads)
#[1] 1772
#[1] "Remaining Features: 7593"
#[1] "Features Removed: 1772"


filtered_dros_head_10percent <- remove.features.10percentcuttoff(proteinGroups_dros_heads)
#[1] 460
#[1] "Remaining Features: 8905"
#[1] "Features Removed: 460"


# mFuzz Soft Clustering----------------------------------------------------

require(Biobase)
exprValues <- new("ExpressionSet", exprs = as.matrix(z))

exprValues.r <- filter.NA(exprValues, thres = 0.9)

exprValues.f = fill.NA(exprValues.r, mode="median")

tmp = filter.std(exprValues.f, min.std = 0.1)

exprValues.s = standardise(tmp)

m1 <- mestimate(exprValues.s)

error <- NA


for(i in 2:50){
  c1 <- mfuzz(exprValues.s, c=i, m=m1)
  error <- rbind(error, c(i,c1$withinerror))
}
plot(error[,1], error[,2])

c1 <- mfuzz(exprValues.s, c=15, m=m1)

pdf(file = "20181022_Phosphosites_SoftClusters.pdf")
par(mar=c(3,3,1,1))
mfuzz.plot2(exprValues.s, cl=c1, mfrow = c(3,2), col.lab="black", x11=F)
dev.off()


imputedHeadGroups <- read_delim ("Timepoints(Branch_ Heads).txt",
                             "\t", escape_double = FALSE, trim_ws = TRUE)

imputedHemolymphGroups <- read_delim ("Timepoint(Branch_ Hemolymph).txt",
                                      "\t", escape_double = FALSE, trim_ws = TRUE)

#count proportion of data that equals zero by column
zero <- colSums(proteinGroups == 0)
zero_head <- zero[17:50]
zero_hemo <- zero[51:58]

#order the columns
zero_head <- zero_head[order(factor(names(zero_head),level = c('LFQ intensity Heads_1', 'LFQ intensity Heads_2','LFQ intensity Heads_3',
                                        'LFQ intensity Heads_4','LFQ intensity Heads_5','LFQ intensity Heads_6',
                                        'LFQ intensity Heads_7','LFQ intensity Heads_8','LFQ intensity Heads_9',
                                        'LFQ intensity Heads_10','LFQ intensity Heads_11','LFQ intensity Heads_12',
                                        'LFQ intensity Heads_13','LFQ intensity Heads_14','LFQ intensity Heads_15',
                                        'LFQ intensity Heads_16','LFQ intensity Heads_17','LFQ intensity Heads_18',
                                        'LFQ intensity Heads_19','LFQ intensity Heads_20','LFQ intensity Heads_21',
                                        'LFQ intensity Heads_22','LFQ intensity Heads_23','LFQ intensity Heads_24',
                                        'LFQ intensity Heads_25','LFQ intensity Heads_26','LFQ intensity Heads_27',
                                        'LFQ intensity Heads_28','LFQ intensity Heads_29','LFQ intensity Heads_30',
                                        'LFQ intensity Heads_31','LFQ intensity Heads_32','LFQ intensity Heads_33',
                                        'LFQ intensity Heads_34')))]


sample = names(zero_hemo)
zerovalues = zero

#create a dataframe with 
dfzero <- data.frame(sample = names(zero_hemo), zerovalues = zero_hemo, totalproteins = rep(10954, each = 8) )
dfzero$additionalproteins <- (dfzero$totalproteins - dfzero$zerovalues)
#bar plots
p <- ggplot(data = dfzero, aes (x = sample, y = zerovalues)) +
        geom_bar(stat = "identity", position = position_dodge())
q <- ggplot(data = dfzero, aes (x = sample, y = totalproteins)) +
  geom_bar(stat = "identity", position = position_dodge())
q + ylim (0, 11000)
p + ylim(0, 11000)

#density plot code
names(proteinGroups)[9] <- "SequenceCoverage" #rename column 9
names(proteinGroups)[1] <- "ProteinID"#rename column 1



Density_sequenceCoverage <- ggplot(proteinGroups, aes(x = SequenceCoverage), fill = "#0073C2FF") 

Density_sequenceCoverage + geom_density() + geom_vline(aes(xintercept = mean(SequenceCoverage)), linetype = "dashed", size = 0.6)


#Imputation choices
#log2protein = imputedHemolymphGroups[,c(7:14)]
#log2protein_heads = imputedHeadGroups[,c(7:40)]

protein_heads = data.frame(proteinGroups[,c(17:50)]) 
protein_heads <- protein_heads[-which(rowSums(protein_heads) == 0),]
protein_heads[protein_heads == 0] <- NA
sum(is.na(protein_heads)) #55185

log2protein_heads = log2(protein_heads) #transform? in Log2 and transpose


str(log2protein_heads)
remove.features(log2protein_heads)


pc<- pca(protein_heads, nPcs = 10, method = "bpca") #no log2 tranformations
pc_log2 <- pca(log2protein_heads, nPcs = 10, method = "bpca") #log2 transformed prior to PCA function
t_pc_log2 <- pca(t(log2protein_heads), nPcs = 10, method = "bpca")
imputed <- completeObs(pc)
imputed_log2 <- completeObs(pc_log2)


write.table(imputed, "Rmodified_imputedProtein_heads.txt", sep="\t")

#log2protein_heads_imputed <- log2(abs(imputed))

plot(pc_log2)
plot(pc)
plot(loadings(pc_log2))
plot(scores(pc))

slplot(pc_log2, scoresLoadings = c(T,T))

log2protein = as.matrix(log2(log2protein_heads))
log2protein = t(log2protein)
#colnames(log2protein) <- seq(1:7108)
colnames(log2protein) <- seq(1:7647)

#log2 <- proteinGroups[proteinGroups == 0] <- NA
#log2 <- log2protein[-(which(row(is.na(log2)) < 0.5)),]
#miss = which(is.na(log2)) # 13609 missing features
#print(paste(length(miss), "missing points."))
#pc.data = pca(log2, method = "bpca", nPcs = 7)
#data.compl = completeObs(pc.data)

pc.data = pca(log2protein, nPcs = 20)

plot(scores(pc.data), pch = 16, main = "PCA Scores Proteins 20181026")
text(scores(pc.data)[,1], scores(pc.data)[,2], labels = rownames(log2protein))

     
  ggplot(proteinGroups, aes(x=SequenceCoverage)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")

ggplot(proteinGroups, aes(x=ProteinID, y=SequenceCoverage) ) + geom_bin2d() +theme_bw()

ggplot(proteinGroups, aes(x=ProteinID, y=SequenceCoverage) ) + geom_density_2d()

