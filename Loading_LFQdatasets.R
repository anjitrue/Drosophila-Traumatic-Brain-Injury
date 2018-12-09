##################################Load Packages #################################################
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
install.packages('GOplot')
devtools::install_github("xia-lab/MetaboAnalystR", build_vignettes=TRUE)
BiocManager::install("ComplexHeatmap", version = "3.8")
BiocManager::install("topGO", version = "3.8")
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
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(topGO)
library(GOplot)
library(cluster)

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

filter.std.plotting <- function (eset, min.std,visu=TRUE){
  #index <- logical(dim(exprs(eset))[1])
  
  tmp <- logical(dim(exprs(eset))[1])
  if (is.numeric(min.std))
  { 
    data <- exprs(eset)
    for (i in 1:length(tmp))
    {
      tmp[i]   <- sd(data[i,],na.rm=TRUE)
      #index[i]  <- ( tmp[i] > min.std)
      
    }
    index <- tmp > min.std
    index[is.na(index)] <- TRUE
    cat(paste(sum(!index),"Proteins excluded.\n"))
  }
  
  if (visu)
  {
    plot(sort(tmp),xlab="Ordered Hemo Proteins",ylab="Standard Deviation")
  }
  eset[index,]
}

fuzzyprep_imputation_included <- function(z){
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  # exclude proteins that have more than 50% of measurements missing
  exprValues.r <- filter.NA(exprValues, thres = 0.50)
  # Fuzzy c-means does not allow for missing values, replace missing values by median values
  exprValues.f = fill.NA(exprValues.r, mode="median")
  #
  tmp = filter.std(exprValues.f, min.std = 0.1)
  # Clustering is performed in Eculidian space, standaridize abundance values to have a mean value of zero
  # Ensures that proteins with similar changes in abundance are close in Euclidean space
  exprValues.s = standardise(tmp)
  return(exprValues.s)
}

fuzzyprep_usepreviousImputation <- function(z){
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  tmp = filter.std.plotting(exprValues, min.std = 0.1, visu = TRUE)
  exprValues.s = standardise(exprValues)
  return(exprValues.s)
}

fuzzyprep_usepreviousImputation_foldchange <- function(z)
{
  exprValues.s <- new("ExpressionSet", exprs = as.matrix(z))
  tmp = filter.std.plotting(exprValues.s, min.std = 0.1)
  #exprValues.s = standardise(exprValues)
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



# set up for baysian pca imputation
hemo_50percent <- as.matrix(filtered_dros_hemo_50percent[,c(4:11)])
hemo_50percent_log2 <- log2(hemo_50percent)
rownames(hemo_50percent_log2) <- filtered_dros_hemo_50percent$id

pc_hemolog2 <- pca(t(hemo_50percent_log2), nPcs = 3, method = "bpca") #pca method
imputed_hemo <- completeObs(pc_hemolog2)

palette(c("mediumorchid2","mediumturquoise","olivedrab3", "darkgoldenrod1", 
          "hotpink3", "red2", "steelblue2", "sienna2","slategray4", 
          "deepskyblue", "orangered", "midnightblue"))

plot(scores(pc_hemolog2), pch = 16,  main = "PCAScores Drosophila Hemolympn")
text(scores(pc_hemolog2)[,1], scores(pc_hemolog2)[,2], labels = colnames(hemo_50percent_log2))

# split data up between the controls and the TBI samples
even_index <- seq(2,8,2) #TBI samples
odd_index <- seq(1,8,2) #controls
imputed_hemo_even <- imputed_hemo[,even_index]
imputed_hemo_odd <- imputed_hemo[,odd_index]
imputed_hemo_reorder <- imputed_hemo[,c(odd_index,even_index)]

#Fold change
imputed_hemo_fold <- imputed_hemo_reorder[,5:8]-imputed_hemo_reorder[,1:4]


# subset for other types of imputation
z_hemo <- proteinGroups_dros_hemo[,4:11]
z_hemo <- log2(z_hemo)
rownames(z_hemo) <- proteinGroups_dros_hemo$id
Z_hemo_even <- z_hemo[,even_index]
z_hemo_odd <- z_hemo[,odd_index]

# subset features with fold change greater than...
important_features_hemo_foldchange <- imputed_hemo_fold[rowSums(abs(imputed_hemo_fold)>1.2)>0,]

# Hemo mFuzz Soft Clustering----------------------------------------------------
z = important_features_hemo_foldchange # change accordingly to use in mfuzzy setup funtions
z = imputed_hemo_reorder
exprValues.s_imput <- fuzzyprep_imputation_included(Z)
exprValues.s_bayesian <- fuzzyprep_usepreviousImputation(z)
exprValues.s_bayesian_foldchange <- fuzzyprep_usepreviousImputation_foldchange(z)
#exprValues.s.nominstd = standardise(exprValues.r) #use data that is non std



error <- NA

y = exprValues.s_bayesian#change according to version of data to perform fuzzy clustering
m1 <- mestimate(y) 


for(i in 2:15){
  c1 <- mfuzz(y, c=i, m=m1)
  error <- rbind(error, c(i,c1$withinerror))
}

plot(error[,1], error[,2])

set.seed(123)
c1 <- mfuzz(y, c=6, m=m1)
c2 <- mfuzz(y, c=4, m=m1)
#pdf(file = "20181022_Phosphosites_SoftClusters.pdf")
#par(mar=c(3,3,1,1))
mfuzz.plot2(y, cl=c1, mfrow = c(3,2), time.labels = c("Ctr_1hr","Ctr_4hr", "Ctr_8hr", "Ctr_24hr","1hr","4hr", "8hr", "24hr"), 
            ylab = "Protein Changes", xlab = "Hemo Time Points",  col.lab="black", x11=F)
mfuzz.plot2(y, cl=c1, mfrow = c(3,2), min.mem = 0.70, time.labels = c("1hr","4hr", "8hr", "24hr","1hr","4hr", "8hr", "24hr"), 
            ylab = "Protein Changes", xlab = "Hemo Time Points",  col.lab="black", x11=F)

mfuzz.plot2(y, cl=c2, mfrow = c(4,2), time.labels = c("1hr","4hr", "8hr", "24hr"), 
            ylab = "Protein Changes", xlab = "Hemo Time Points",  col.lab="black", x11=F)
mfuzz.plot2(y, cl=c2, mfrow = c(4,2), time.labels = c("1hr","4hr", "8hr", "24hr"), 
            min.mem = 0.50, ylim.set = c(-2,2), ylab = "Protein Changes", xlab = "Hemo Time Points",  col.lab="black", x11=F)

### Fuzzy membership filter
#compare two feature clusters (fold change*)

table(rowSums(c1$membership > .70) >0)

#dev.off()

a_bayesian <- acore(exprValues.s_bayesian_foldchange, c2, min.acore=0.0)

cluster1 <- data.frame(a_bayesian[[1]])
cluster2 <- data.frame(a_bayesian[[2]])
cluster3 <- data.frame(a_bayesian[[3]])
cluster4 <- data.frame(a_bayesian[[4]])
#For Phospho Data
cluster5 <- data.frame(a_bayesian[[5]])
cluster6 <- data.frame(a_bayesian[[6]])



names.FromCluster_2 <- cluster2$NAME
cluster3_imputed_hemo_reorder <- imputed_hemo_reorder[rownames(imputed_hemo_reorder) %in% names.FromCluster, ]
imputed_hemo_cluster3_fold <- as.matrix(cluster3_imputed_hemo_reorder[,5:8]-cluster3_imputed_hemo_reorder[,1:4])

rownames(proteinGroups_dros_hemo) <- proteinGroups_dros_hemo$id
cluster2_hemo_proteins <-  proteinGroups_dros_hemo[rownames(proteinGroups_dros_hemo) %in% names.FromCluster,]
genes_cluster2 <- cluster2_hemo_proteins$Gene.names
write.csv(cluster2_hemo_proteins, "cluster2_hemo_proteins.csv")

#distance clustering reordered hemolympth data
#dis <- dist(imputed_hemo_cluster2_fold,method = "euclidean")
#distance = dist(imputed_hemo_cluster2_fold, method = "euclidean")
#row.cluster = hclust(imputed_hemo_cluster2_fold, method = "aver")


reorder_hclust <- hclust(dist(imputed_hemo_reorder[rowSums(c1$membership >0.7)>0,], method = "euclidean"))

clusters <- cutree(reorder_hclust, k = 6)

par(mar=c(3,4,3,3)+0.1)
plot(reorder_hclust, label= FALSE)

rect.hclust(reorder_hclust, k=6, border = "red")

scaleRYG <- colorRampPalette(c("red","black","darkgreen"), space = "rgb")(31)
breaksList = c(-10,seq(-2, 2, length.out = 28),10)

pheatmap(
  mat               = imputed_hemo_reorder[rowSums(c1$membership >0.7)>0,],
  color             = scaleRYG,
  border_color      = NA,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 10,
  main              = "Default Heatmap"
)

table(rowSums(abs(imputed_hemo_fold)>1.2)>0)
p <- pheatmap(
  mat               = imputed_hemo_fold[rowSums(abs(imputed_hemo_fold)>1.2)>0,], #[rowSums(c2$membership >0.5)>0,],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = c("1hr","4hr", "8hr", "24hr"),
  breaks = breaksList,
  clustering_distance_rows = "euclidean",
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 10,
  main              = "Fold Change relative to Control - Hemo"
)
dev.off()

important_features_foldchange_compress <- important_features_hemo_foldchange
important_features_foldchange_compress[important_features_foldchange_compress < -2] = -2.1
important_features_foldchange_compress[important_features_foldchange_compress > 2] = 2.1


pheatmap(
  mat               = important_features_foldchange_compress, #[rowSums(c2$membership >0.5)>0,],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = c("1hr","4hr", "8hr", "24hr"),
  #breaks = breaksList,
  clustering_distance_rows = "euclidean",
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 10,
  main              = "Fold Change with compressed Hemo Data, binning <-2 and >2"
)
dev.off()
fold_hclust <- hclust(dist(important_features_foldchange_compress, method = "euclidean"))
fold_hclust <- agnes(dist(imputed_hemo_fold[rowSums(abs(imputed_hemo_fold)>1.2)>0,], method = "euclidean"))

clusters <- cutree(fold_hclust, k = 4) #k = 20)

plot(fold_hclust, label= FALSE)

#rect.hclust(fold_hclust, k = 20, border = "red")
rect.hclust(fold_hclust, k=3, border = "red")


cluster2 <- kmeans(t(as.matrix(imputed_hemo_fold[rowSums(abs(imputed_hemo_fold)>1.2)>0,])),10)

clusplot(important_features_foldchange_compress, clusters, lines = 0)
clusplot(imputed_hemo_fold[rowSums(abs(imputed_hemo_fold)>1.2)>0,], cluster2, lines = 0)
clusplot(as.matrix(reorder_hclust), clusters, lines = 0)

write.csv(cluster1, "HemoOrdered_bayesian_cluster1Info.csv")
write.csv(cluster2, "HemoOrdered_bayesian_cluster2Info.csv")
write.csv(cluster3, "HemoOrdered_bayesian_cluster3Info.csv")
write.csv(cluster4, "HemoOrdered_bayesian_cluster4Info.csv")
write.csv(cluster5, "HemoOrdered_bayesian_cluster5Info.csv")
write.csv(cluster6, "HemoOrdered_bayesian_cluster6Info.csv")



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

