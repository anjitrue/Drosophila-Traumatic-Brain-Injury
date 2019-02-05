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
library(yaml)
library(factoextra)
library(ellipse)


##### Load data #####
proteinGroups_dros_heads <- read.csv("E:/Projects/Proteomics/DorsophilaHead_Experiment/txt_dros_heads_plusFrac/proteinGroups__dros_heads_plusFrac.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)
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

#log2 transform
LFQheads_LOG2_complete <- log2(as.matrix(LFQheads_complete[,c(4:37)]))


#set row names to match protein group identifer number
rownames(LFQheads_LOG2_complete) <- LFQheads_complete$id

samples <- as.numeric(sub(".*Heads_","", colnames(LFQheads_LOG2_complete)))
colnames(LFQheads_LOG2_complete) <- samples
LFQheads_LOG2 <- LFQheads_LOG2_complete[,order(as.numeric(as.character(colnames(LFQheads_LOG2_complete))))]

#rename the column names that have been sorted numerically
even = seq(2,34,2)
odd = seq(1,33,2)
samples <- seq(1:34)
colnames(LFQheads_LOG2) <- samples


head_meta <- data.frame(Samples = samples,Sample_Type=rep(c("Control","TBI"),17))
head_meta$time <- c(rep("day", 20), rep("2 days",2), rep("3 days",2), rep("weeks", 10))

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

colors <- c("#F47F72", #coral
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
sample.colors = as.numeric(factor(head_meta$time))

plot(pca_complete_LFQheads$x, pch = 19, col = sample.colors, main = "Principle Component Analysis\n Complete Head Data", ylim = c(-50, 50), xlim = c(-50, 50), xlab = "PC1\n32.48%", ylab ="PC2\n10.27%")
text(pca_complete_LFQheads$x[,1], pca_complete_LFQheads$x[,2], labels = colnames(LFQheads_LOG2), 
     col = sample.colors)
legend("topright", legend = levels(factor(head_meta$time)), pch = 16, col = 1:length(levels(factor(head_meta$time))), y.intersp = 0.7)


# plot PCA with samples within the first 24 hours
sample.colors = as.numeric(factor(head_meta[1:20,]$Sample_Type))
plot(pca_first24hours$x, pch = 19, cex = 2, col = sample.colors, main = "Principle Component Analysis\n First 24 hours", ylim = c(-50, 50), xlim = c(-30, 30), xlab = "PC1\n 17.20%", ylab ="PC2\n 14.42%")
text(pca_first24hours$x[,1], pca_first24hours$x[,2], labels = colnames(LFQheads_LOG2[,1:20]))
legend("topleft", legend = levels(factor(head_meta[1:20,]$Sample_Type)), pch = 16, col = 1:length(levels(factor(head_meta[1:20,]$Sample_Type))))

# plot PCA with samples across multiple days
sample.colors = as.numeric(factor(head_meta[21:34,]$Sample_Type))
plot(pca_days$x, pch = 19, cex = 2, col = sample.colors, main = "Principle Component Analysis\n Across Multiple Days", ylim = c(-50, 50), xlim = c(-30, 30), xlab = "PC1\n 16.31%", ylab ="PC2\n 9.57%")
text(pca_days$x[,1], pca_days$x[,2], labels = colnames(LFQheads_LOG2[,21:34]))
legend("topleft", legend = levels(factor(head_meta[21:34,]$Sample_Type)), pch = 16, col = 1:length(levels(factor(head_meta[21:34,]$Sample_Type))))



###### Plot of variance explained by PCA ######

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot1 <- barplot(summary(pca_complete_LFQheads)$importance[2,1:34]*100, ylab ="percent variance explained", main = "PCA scree plot", ylim=c(0,100), col = colors[6])
lines(x = plot1, y =summary(pca_complete_LFQheads)$importance[3,1:34]*100, col = colors[5], pch=19, type = "b")
legend("topleft", pch=19, col=colors[5], "cummulative variance", bty="n")

##### PCA_alltimes with confidence circles
par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
#define the levels associated with your samples, in this case Control and TBI 
lvs <- levels(factor(head_meta$time))
#create an array that is an arbritrary length of 100 by 2, and 3rd length is defined as the number of arrays for the defined levels of the data
pts.array <- array(0, dim = c(100, 2, length(lvs)))
#for loop, iterate through the levels of your data
i= 3
for (i in 1:length(lvs)) 
{
  # create a variable that determines which samples are in the current level 
  inx <- head_meta$time == lvs[i]
  # calculate variation of your data in current level
  groupVar <- var(cbind(pca_complete_LFQheads$x[inx,1], pca_complete_LFQheads$x[inx,2]), na.rm = T)
  # calculate group mean of your data in current level
  groupMean <- cbind(mean(pca_complete_LFQheads$x[inx,1], na.rm = T), mean(pca_complete_LFQheads$x[inx,2], na.rm = T))
  
  #fill in array for the current level
  pts.array[, , i] <- ellipse(groupVar, centre = groupMean, 
                              level = 0.95, npoints = 100)
}

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
plot(pca_complete_LFQheads$x, col= colors[(2:6)][as.factor(head_meta$time)], xlab = "PC1\n32.48", ylab ="PC2\n 10.27", pch=19, cex =2, main= "Principle Component Analysis", ylim = c(-50, 50), xlim = c(-50, 50))
legend("topleft", col = colors[2:6], pch=19, bty = "n", legend =levels(as.factor(head_meta$time)))

# polygon(pts.array[, , 1], col = adjustcolor(colors[1], alpha = 0.25), border = NA)
# 
# polygon(pts.array[, , 2], col = adjustcolor(colors[2], alpha = 0.25), border = NA)
# 
# polygon(pts.array[, , 3], col = adjustcolor(colors[3], alpha = 0.25), border = NA)
# 
# polygon(pts.array[, , 4], col = adjustcolor(sample.colors[4], alpha = 0.25), border = NA)


text(pca_complete_LFQheads$x[,1], pca_complete_LFQheads$x[,2], samples, pos =4, offset = 0.5, cex = 0.8)


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
  mat               = important_fold_change_brains[,11:17],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  cutree_row = 7,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = times[11:17],
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
  mat               = compress_fold_brains [,11:17],
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  cutree_row = 7,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = times[11:17],
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

