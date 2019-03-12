install.packages("BiocManager")
BiocManager::install("pcaMethods", version = "3.8")
BiocManager::install("Mfuzz", version = "3.8")

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
library(cluster)
library(pls)
library(ellipse)


##### Load data #####
proteinGroups_dros_heads <- read.csv("G:/Projects/Proteomics/DorsophilaHead_Experiment/txt_dros_heads_plusFrac/proteinGroups__dros_heads_plusFrac.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
cat(paste0("Number of protein groups in Heads Data before any filtering of missing measurements: ", nrow(proteinGroups_dros_heads))) #9381


##### Format data #####

subsetLFQ <- function(q){
  y <- q[,which(names(q) %in% c("Protein.IDs", "id","Gene.names"))]
  z <- q[,grep("LFQ.intensity",names(q))]
  z[z == 0] <- NA
  x <- bind_cols(y,z)
  #how many missing values in protein groups dataframe
  print(paste("Number of missing measurements", sum(is.na(x)), " out of ", ncol(x)*nrow(x)))
  print(paste("% of missing data", (sum(is.na(x)))/(ncol(x)*nrow(x))*100))
  x <- x[complete.cases(x),]
  return(x)
}

#Subset LFQ data to complete cases only
LFQheads_complete <- subsetLFQ(proteinGroups_dros_heads)
dim(LFQheads_complete) #[1] 5507   37

# # TO DO plot data of what is missing
# remove.features.50percentcuttoff <- function (x) {
#   features.missing = rowMeans(is.na(x)) 
#   print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
#   features.missing.50more = rownames(x)[features.missing > 0.50] 
#   
#   keep.features = which(features.missing <= 0.50) 
#   print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
#   names(keep.features) = keep.features 
#   
#   remove.features = which(features.missing > 0.50)
#   print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
#   names(remove.features) = remove.features
#   
#   filtered = x[-which(rownames(x) %in% remove.features),]
#   return(filtered)
# }

#log2 transform and subset only complete cases
LFQheads_LOG2_complete <- log2(as.matrix(LFQheads_complete[,c(4:37)]))


#set row names to match protein group identifer number
rownames(LFQheads_LOG2_complete) <- LFQheads_complete$id

samples <- as.numeric(sub(".*Heads_","", colnames(LFQheads_LOG2_complete)))
colnames(LFQheads_LOG2_complete) <- samples
LFQheads_LOG2 <- LFQheads_LOG2_complete[,order(as.numeric(as.character(colnames(LFQheads_LOG2_complete))))]

#rename the column names that have been sorted numerically
even = seq(2,34,2)
odd = seq(1,33,2)
samples <- seq(1:34) #currently colnames are characters, switch to numeric
colnames(LFQheads_LOG2) <- samples


head_meta <- data.frame(Samples = samples,Sample_Type=rep(c("Control","TBI"),17))
head_meta$day <- c(rep("day", 20), rep("2 days",2), rep("3 days",2), rep("weeks", 10))
number_times = c(0, 0, (.5/24), (.5/24), (1/24), (1/24), (2/24), (2/24), 
                 (4/24), (4/24), (6/24), (6/24), (8/24), (8/24), (12/24), (12/24), 
                 (16/24), (16/24), (24/24), (24/24), (48/24), (48/24), (72/24), (72/24), 
                 7, 7, 14, 14, 21, 21, 28, 28, 35, 35)
head_meta$time <- number_times
##### Plotting###

#### PCA prcomp function complete cases ####
pca_complete_LFQheads <- prcomp(t(LFQheads_LOG2))

pca_first24hours <- prcomp(t(LFQheads_LOG2[,1:20]))

pca_days <- prcomp(t(LFQheads_LOG2[,21:34]))

# Print standard deviations of the PC's determined
print(pca_complete_LFQheads)


# Plot function outputs a plot with the variance of the data on the y-axis and PC on x-axis
plot(pca_complete_LFQheads, type = "l")

#Find percentages here
summary(pca_complete_LFQheads)

summary(pca_first24hours)

summary(pca_days)

colors <- c("#206F94", # teal
            "#F47F72", #coral
            "#75C69D", #baby green
            "#2CA8E0", #coon blue
            "#1F6F94", #darker blue
            "#5C66AF", #light purpleblue
            "#2A4093", #dark purpleblue
            "#2C9379", #dark baby green
            "#83D5F7", #light coon blue
            "#93211E", #dark red
            "#E73C25", #red
            "#81143A", #dark pink
            "#ED237A") #hot pink)

# plot PCA with all samples coloring by time of collection
sample.colors = as.numeric(factor(head_meta$day))
plot(pca_complete_LFQheads$x, pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Complete Head Data", 
     ylim = c(-50, 50), xlim = c(-50, 50), 
     xlab = "PC1\n32.48%", ylab ="PC2\n10.27%")
text(pca_complete_LFQheads$x[,1], pca_complete_LFQheads$x[,2], 
     labels = colnames(LFQheads_LOG2), 
     col = colors[sample.colors])
legend("topright", legend = levels(factor(head_meta$day)), pch = 16, 
       col = colors[1:length(levels(factor(head_meta$day)))], y.intersp = 0.7)


# plot PCA with samples within the first 24 hours
sample.colors = as.numeric(factor(head_meta[1:20,]$Sample_Type))
plot(pca_first24hours$x, pch = 19, cex = 3, col = colors[sample.colors], 
     main = "Principle Component Analysis\n First 24 hours", 
     ylim = c(-25, 25), xlim = c(-25, 25), 
     xlab = "PC1\n 17.20%", ylab ="PC2\n 14.42%")
text(pca_first24hours$x[,1], pca_first24hours$x[,2], 
     labels = colnames(LFQheads_LOG2[,1:20]))
legend("topleft", legend = levels(factor(head_meta[1:20,]$Sample_Type)), 
       pch = 16, 
       col = colors[1:length(levels(factor(head_meta[1:20,]$Sample_Type)))])

# PC3, PC2
sample.colors = as.numeric(factor(head_meta[1:20,]$Sample_Type))
plot(pca_first24hours$x[,3], pca_first24hours$x[,2], pch = 19, cex = 3, 
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n First 24 hours", 
     ylim = c(-25, 25), xlim = c(-25, 25), 
     xlab = "PC3 9.07%", ylab ="PC2 14.42%")
text(pca_first24hours$x[,3], pca_first24hours$x[,2], 
     labels = colnames(LFQheads_LOG2[,1:20]))
legend("topleft", legend = levels(factor(head_meta[1:20,]$Sample_Type)), 
       pch = 16, 
       col = colors[1:length(levels(factor(head_meta[1:20,]$Sample_Type)))])

# plot PCA with samples across multiple days
sample.colors = as.numeric(factor(head_meta[21:34,]$Sample_Type))
plot(pca_days$x[,3],pca_days$x[,2], pch = 19, cex = 3, 
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Across Multiple Days", 
     ylim = c(-25, 25), xlim = c(-25, 25), 
     xlab = "PC3 10.16%", ylab ="PC2 11.49%")
text(pca_days$x[,3], pca_days$x[,2], labels = colnames(LFQheads_LOG2[,21:34]))
legend("topleft", legend = levels(factor(head_meta[21:34,]$Sample_Type)), 
       pch = 16, 
       col = colors[1:length(levels(factor(head_meta[21:34,]$Sample_Type)))])



###### Plot of variance explained by PCA ######
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot1 <- barplot(summary(pca_complete_LFQheads)$importance[2,1:34]*100, ylab ="percent variance explained", main = "PCA scree plot", ylim=c(0,100), col = colors[6])
lines(x = plot1, y =summary(pca_complete_LFQheads)$importance[3,1:34]*100, col = colors[5], pch=19, type = "b")
legend("topleft", pch=19, col=colors[5], "cummulative variance", bty="n")

##### PCA_all time points with confidence circles #####

#define the levels associated with your samples, in this case Control and TBI 
lvs <- levels(factor(head_meta$day))
#create an array that is an arbritrary length of 100 by 2, and 3rd length is defined as the number of arrays for the defined levels of the data
pts.array <- array(0, dim = c(100, 2, length(lvs)))
#for loop, iterate through the levels of your data
i = 3
for (i in 1:length(lvs)) 
{
  # create a variable that determines which samples are in the current level 
  inx <- head_meta$day == lvs[i]
  # calculate variation of your data in current level
  groupVar <- var(cbind(pca_complete_LFQheads$x[inx,1], pca_complete_LFQheads$x[inx,2]), na.rm = T)
  # calculate group mean of your data in current level
  groupMean <- cbind(mean(pca_complete_LFQheads$x[inx,1], na.rm = T), mean(pca_complete_LFQheads$x[inx,2], na.rm = T))
  
  #fill in array for the current level
  pts.array[, , i] <- ellipse(groupVar, centre = groupMean, 
                              level = 0.95, npoints = 100)
}

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_complete_LFQheads$x, col= colors[(2:6)][as.factor(head_meta$day)], xlab = "PC1\n32.48", ylab ="PC2\n 10.27", pch=19, cex =2, main= "Principle Component Analysis", ylim = c(-50, 50), xlim = c(-50, 50))
legend("topleft", col = colors[2:6], pch=19, bty = "n", legend =levels(as.factor(head_meta$day)))

polygon(pts.array[, , 1], col = adjustcolor(colors[1], alpha = 0.25), border = NA)
# 
polygon(pts.array[, , 2], col = adjustcolor(colors[2], alpha = 0.25), border = NA)
# 
polygon(pts.array[, , 3], col = adjustcolor(colors[3], alpha = 0.25), border = NA)
# 
polygon(pts.array[, , 4], col = adjustcolor(sample.colors[4], alpha = 0.25), border = NA)


text(pca_complete_LFQheads$x[,1], pca_complete_LFQheads$x[,2], samples, pos =4, offset = 0.5, cex = 0.8)


##### PCA for Days After samples ######
lvs <- levels(factor(head_meta$Sample_Type))
pts.array <- array(0, dim = c(100, 2, length(lvs)))

for(i in 1:length(lvs)){
  inx <- head_meta[21:34,]$Sample_Type == lvs[i]
  groupVar <- var(cbind(pca_days$x[inx,3], pca_days$x[inx,2]), na.rm = T)
  groupMean <- cbind(mean(pca_days$x[inx,3], na.rm = T), mean(pca_days$x[inx,2], na.rm = T))
  pts.array[,,i] <- ellipse(groupVar, centre = groupMean, level = 0.95, npoints = 100)
}

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_days$x[,3], pca_days$x[,2], col = colors[1:2][as.factor(head_meta[21:34,]$Sample_Type)], 
     xlab = , ylab = , pch=19, cex=2, main = "Principle Component Analysis", ylim = c(-40,40), xlim = c(-30,30))
polygon(pts.array[, , 1], col = adjustcolor(colors[1], alpha = 0.25), border = NA)
polygon(pts.array[, , 2], col = adjustcolor(colors[2], alpha = 0.25), border = NA)


###### PCA Loadings #####
#All samples
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_complete_LFQheads$rotation, xlim = c(-.10,.20),
     xlab = "PC1 Loadings" , ylab = "PC2 Loadings", col = "#5C66AF",  
     pch=19, cex=2, main = "Loadings from Principle Component Analysis")
text(pca_complete_LFQheads$rotation[,1], pca_complete_LFQheads$rotation[,2], rownames(LFQheads_LOG2), pos = 4, offset = 0.5, cex = 0.8)

#First 24 hours PC1 vs PC2
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_first24hours$rotation, xlim = c(-.26,.15),
     xlab = "PC1 Loadings" , ylab = "PC2 Loadings", col = "#5C66AF",  
     pch=19, cex=2, main = "Loadings from Principle Component Analysis")
text(pca_first24hours$rotation[,1], pca_first24hours$rotation[,2], rownames(LFQheads_LOG2), pos = 4, offset = 0.5, cex = 0.8)

#PC3 vs PC2
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_first24hours$rotation[,3], pca_first24hours$rotation[,2], xlim = c(-.26,.15),
     xlab = "PC3 Loadings" , ylab = "PC2 Loadings", col = "#5C66AF",  
     pch=19, cex=2, main = "Loadings from Principle Component Analysis")
text(pca_first24hours$rotation[,3], pca_first24hours$rotation[,2], rownames(LFQheads_LOG2), pos = 4, offset = 0.5, cex = 0.8)

#Days after samples
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_days$rotation, xlim = c(-.15,.20),
     xlab = "PC1 Loadings" , ylab = "PC2 Loadings", col = "#5C66AF",  
     pch=19, cex=2, main = "Loadings from Principle Component Analysis")
text(pca_days$rotation[,1], pca_days$rotation[,2], rownames(LFQheads_LOG2), pos = 4, offset = 0.5, cex = 0.8)

#PC3 vs PC2
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_days$rotation[,3], pca_days$rotation[,2], 
     xlab = "PC3 Loadings" , ylab = "PC2 Loadings", col = "#5C66AF",  pch=19, cex=2, main = "Loadings from Principle Component Analysis")
text(pca_days$rotation[,3], pca_days$rotation[,2], rownames(LFQheads_LOG2), pos = 4, offset = 0.5, cex = 0.8)

##### Pls-da of brain #####

#By treatement

treatment_plsda <- plsr(as.numeric(head_meta$Sample_Type) ~ t(LFQheads_LOG2), method = "oscorespls", ncomp =  4)
summary(treatment_plsda)

sample.colors = as.numeric(factor(head_meta$Sample_Type))
plot(treatment_plsda$scores, pch = 19,  col = colors[sample.colors])
legend("topleft", legend = levels(factor(head_meta$Sample_Type)), pch = 16, col = colors[1:length(levels(factor(head_meta$Sample_Type)))])

lvs <- levels(head_meta$Sample_Type)
pts.array <- array(0, dim = c(100,2, length(lvs)))
for(i in 1:length(lvs)){
  inx <- head_meta$Sample_Type == lvs[i]
  groupVar <- var(cbind(treatment_plsda$scores[inx,1], treatment_plsda$scores[inx,2]), na.rm = T)
  groupMean <- cbind(mean(treatment_plsda$scores[inx,1]), mean(treatment_plsda$scores[inx,2]))
  pts.array[,,i] <- ellipse(groupVar, center = groupMean, level = 0.95, npoints = 100)
}

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(treatment_plsda$scores, col = colors[(1:2)][as.factor(head_meta$Sample_Type)], xlab = "Component 1" , ylab = "Component 2",
     pch = 19, main = " PLS-DA", cex = 2)
polygon(pts.array[,,1], col = adjustcolor(colors[1], alpha = 0.25), border = NA)
polygon(pts.array[,,2], col = adjustcolor(colors[2], alpha = 0.25), border = NA)
text(treatment_plsda$scores[,1], treatment_plsda$scores[,2], as.character(head_meta$Samples), pos = 4, offset = 0.5, cex = 0.8)
legend("topleft", legend = levels(factor(head_meta$Sample_Type)), pch = 16, col = colors[1:length(levels(factor(head_meta$Sample_Type)))])

dev.off()

dotchart(rev(treatment_plsda$loadings[order(-abs(treatment_plsda$loadings[,1])),1][1:25]), 
         pch = 19,
         color = "#5C66AF",
         xlab = "Component 1 Loadings")
plot(treatment_plsda$loadings[,1], pch = 19, col = "#5C66AF")
text(treatment_plsda$loadings[,1], as.character(gsub(".*)", "", rownames(treatment_plsda$loadings))))


plot(treatment_plsda$loadings[order(treatment_plsda$loadings[,1]),1])

##### Clustering in HeatMap #####

reorder_brains <- LFQheads_LOG2[,c(odd, even)]
times <- c("0min", "30min", "1hr", "2hr", "4hr", "6hr", "8hr", "12hr", "16hr", "24hr", "48hr", "72hr", "7days", "14days", "21days", "28days", "35days")

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
boxplot(reorder_brains, notch = TRUE)

#histogram of log2(intensities) for data set that has NA's
hist(LFQheads_LOG2, breaks = 20, xlab = "Log2(Brain data complete cases)", main = "Histogram of Brain Data")

scaleRYG <- colorRampPalette(c("red","black","chartreuse4"), space = "rgb")(31)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
pheatmap(
  mat               = reorder_brains,
  color             = scaleRYG,
  border_color      = NA,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  #labels_col = c(sprintf("Ctr_%d", seq(1:17)), sprintf("TBI_%d", seq(1:17))),
  labels_col = sub(" ", "_", c(paste("Ctr",times), paste("TBI", times))),
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "Brains LOG2(Intensity)"
)

##### Fold Change #####

fold_change_brains <- reorder_brains[,18:34] - reorder_brains[,1:17]
summary(fold_change_brains)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
boxplot(fold_change_brains, notch = TRUE)

important_fold_change_brains <- fold_change_brains[rowSums(abs(fold_change_brains)>1.2)>0,]
dim(important_fold_change_brains) #[1] 809  17

compress_fold_brains <- important_fold_change_brains
compress_fold_brains[compress_fold_brains < -4] = -4.1
compress_fold_brains[compress_fold_brains > 4] = 4.1

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
pheatmap(
  mat               = important_fold_change_brains[,1:10],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  #cutree_row = 7,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = times[1:10],
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "Fold Change first Days Later - Brains"
)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
pheatmap(
  mat               = compress_fold_brains [,1:10],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  #cutree_row = 7,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = times[1:10],
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "Compressed FC First Days Later - Brains, binning <-4 and >4"
)

# Proteins with fold change greater than 1.2:
table(rowSums(abs(fold_change_brains)>1.2)>0) 

# FALSE  TRUE 
# 4698   809

# Protein with fold change greater than 2:
table(rowSums(abs(fold_change_brains)>2)>0)

# FALSE  TRUE 
# 5245   262

# Protein with fold change greater than 4:
table(rowSums(abs(fold_change_brains)>4)>0)

# FALSE  TRUE 
# 5467    40

# Protein with fold change greater than 8:
table(rowSums(abs(fold_change_brains)>8)>0)

# FALSE  TRUE 
# 5505     2

##### Eucclidean clustering of Fold Change #####
LFQheads_LOG2_hclust <- hclust(dist(LFQheads_LOG2[,1:20], method = "euclidean"))
clusters <- cutree(LFQheads_LOG2_hclust, k=20)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(LFQheads_LOG2_hclust, label= FALSE, main = "Brain Dendogram Days After - Euclidean Distance")
rect.hclust(LFQheads_LOG2_hclust, k=10, border = "red")

#fold change clustering
FC_brains_hclust <- hclust(dist(fold_change_brains, method = "euclidean"))
clusters <- cutree(FC_brains_hclust, k =20)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(FC_brains_hclust, label= FALSE, main = "Brain Dendogram FC - Euclidean Distance")
rect.hclust(FC_brains_hclust, k=20, border = "red")

filter.std.plotting <- function (eset, min.std,visu=TRUE)
{
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
    cat(paste(sum(!index),"Proteins have a standard deviation greater than ", min.std, ".\n"))
  }
  
  if (visu)
  {
    plot(sort(tmp),xlab="Ordered Hemo Proteins",ylab="Standard Deviation")
  }
  eset[index,]
}

fuzzyprep_imputation_included <- function(z)
{
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  # exclude proteins that have more than 50% of measurements missing
  exprValues.r <- filter.NA(exprValues, thres = 0.50)
  # Fuzzy c-means does not allow for missing values, replace missing values by median values
  exprValues.f = fill.NA(exprValues.r, mode="median")
  # Set a minimum threshold for variation 
  tmp = filter.std(exprValues.f, min.std = 0.1)
  # Clustering is performed in Eculidian space, standaridize abundance values to have a mean value of zero
  # Ensures that proteins with similar changes in abundance are close in Euclidean space
  exprValues.s = standardise(tmp)
  return(exprValues.s)
}

fuzzyprep_usepreviousImputation <- function(z)
{
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  tmp = filter.std.plotting(exprValues, min.std = 0.1)
  exprValues.s = standardise(exprValues)
  return(exprValues.s)
}

fuzzyprep_usepreviousImputation_foldchange <- function(z)
{
  exprValues <- new("ExpressionSet", exprs = as.matrix(z))
  tmp = filter.std.plotting(exprValues, min.std = 0.1)
  #exprValues.s = standardise(exprValues)
  return(exprValues)
}

##### Finding proteins of interest #####
evenTimes <- head_meta[even,]$time
labels_hours <- c("0","0.5","1","2","4","6","8","12","16","24")
labels_days <- c("2","3","7","14","21","28","35")

proteinOfInterest <- function(y){
  print(paste("proteinID: ", LFQheads_complete[which(LFQheads_complete$id == y),1]))
  print(paste("geneID: ", LFQheads_complete[which(LFQheads_complete$id == y),2]))
  plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==y),even], type = "l")
  lines(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==y), odd], col = "red") 
}

# Input is id of protein
proteinOfInterest(n)

###### Time series plots ####

ranked_pls_loading1 <- as.numeric(sub(".*LOG2)","", 
                                      names(rev(treatment_plsda$loadings[order(-abs(treatment_plsda$loadings[,1])),1][1:50]))))

timeSeries_protein_plot <- function(n){
  even_timepoints <- LFQheads_LOG2[which(rownames(LFQheads_LOG2)==n),even]
  odd_timepoints <- LFQheads_LOG2[which(rownames(LFQheads_LOG2)==n),odd]
  
  if(min(even_timepoints) < min(odd_timepoints)){
    minimum_y_Intensity <- min(even_timepoints)
  }else{
    minimum_y_Intensity <- min(odd_timepoints)
  }
  
  if(max(even_timepoints) > max(odd_timepoints)){
    max_y_Intensity <- max(even_timepoints)
  }else{
    max_y_Intensity <- max(odd_timepoints)
  }
  
  #par(mfrow = c(1,2))
  #layout(matrix(c(1,1), 1, 2, byrow = TRUE), widths=c(2,1), heights=c(1,2))
  par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
  #plot for hours
  plot(evenTimes[1:10], even_timepoints[1:10],
       ylim = c(minimum_y_Intensity,max_y_Intensity),
       type = "b,c",
       pch = 19,
       lty = 2, 
       xaxt = "n")
  lines(evenTimes[1:11],odd_timepoints[1:11], col = "red") 
  axis(1,evenTimes[1:10], labels = labels_hours)
  #plot for days
  plot(evenTimes[11:17], even_timepoints[11:17], 
       ylim = c(minimum_y_Intensity,max_y_Intensity), 
       type = "b,c",
       pch = 19,
       lty = 2, 
       xaxt = "n")
  lines(evenTimes[10:17],odd_timepoints[10:17], col = "red")
  axis(1,evenTimes[11:17], labels = labels_days)
}

#plot the time series plot for individual proteins
n <- 3272 #LFQheads_complete[which(LFQheads_complete$id == n),]
timeSeries_protein_plot(n)

#ks.test(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==n),even[10:17]],LFQheads_LOG2[which(rownames(LFQheads_LOG2)==n), odd[10:17]])

df <- data.frame(Protein.ID = as.numeric(), Anova.p.Samply_Type = as.numeric(), 
                 Anova.p.time = as.numeric(), Anova.p.Samply_Type.Time = as.numeric())


for(i in 1:length(ranked_pls_loading1)){
  n <- ranked_pls_loading1[i]
  
  aov_n <- aov(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==n),] ~  head_meta$Sample_Type * head_meta$time)
  
  df[i,] <- c(n,summary(aov_n)[[1]][["Pr(>F)"]])
}



#individual anova calculations 
anova1 <- aov(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==4034),] ~  head_meta$Sample_Type * head_meta$time)
summary(anova1)

anova5606 <- aov(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==5606),] ~  head_meta$Sample_Type * head_meta$time)
summary(anova5606)


anova2535 <- aov(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==2535),] ~  head_meta$Sample_Type * head_meta$time)
summary(anova2535)

#look into K.S. test for best line of fits

help("ks.test")


# 4362 - Loopin-1
LFQheads_complete[which(LFQheads_complete$id == 4362),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==4362),even])

# 5606
LFQheads_complete[which(LFQheads_complete$id == 5606),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==5606),even])


# 2535
LFQheads_complete[which(LFQheads_complete$id == 2535),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==2535),even])

# 8853
LFQheads_complete[which(LFQheads_complete$id == 8853),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==8853),even])

# 3270
LFQheads_complete[which(LFQheads_complete$id == 3270),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==3270),even])


# 1145
x <- LFQheads_complete[which(LFQheads_complete$id == 1145),]
plot(LFQheads_LOG2[which(rownames(LFQheads_LOG2)==1145),even])

# Look for individual proteins by Gene.Name
proteinGroups_dros_heads[which(proteinGroups_dros_heads$Gene.names == "Tau"),]


Protiens <- rownames(proteinGroups_dros_heads[grepl("Q", proteinGroups_dros_heads$Pr),])
df_proteins <- proteinGroups_dros_heads[Protiens,]

x <- rnorm(50)
y <- runif(30)

plot(x,y)
ks.test(x,z)
