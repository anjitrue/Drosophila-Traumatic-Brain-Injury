library(devtools)
library(pcaMethods)
library(ggplot2)
library(reshape2)
library(factoextra)
library(dplyr)
library(pheatmap)

#### Load Data ####

proteinGroups_dros_hemo_MAXQuant <- read.csv("G:/Projects/Proteomics/DorsophilaHead_Experiment/txt_hemo_plusFrac/proteinGroups.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

cat(paste0("Number of protein groups in Hemolymph Data before any filtering of missing measurements: ",
                 nrow(proteinGroups_dros_hemo))) #, "\n", "Number of columns: ", ncol(proteinGroups_dros_hemo)))


#### Format Data ####

subsetLFQ <- function(x){
  y <- proteinGroups_dros_hemo_MAXQuant[,which(names(proteinGroups_dros_hemo_MAXQuant) %in% c("Protein.IDs", "id","Gene.names"))]
  z <- x[,grep("LFQ.intensity",names(x))]
  z[z == 0] <- NA
  x <- data.frame(y,z)
  x <- x[,-grep("Frac",names(x))]
  return(x)
}

remove.features.50percentcuttoff <- function (x) {
  features.missing = rowMeans(is.na(x)) 
  print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  keep.features = which(features.missing <= 0.50) 
  print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
  names(keep.features) = keep.features 
  
  remove.features = which(features.missing > 0.50)
  print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
  names(remove.features) = remove.features
  
  filtered = x[-which(rownames(x) %in% remove.features),]
  return(filtered)
}


#subset proteins to only protein.Id, gene.id, id, and LFQ not including the fractionated sample
proteinGroups_dros_hemo <- subsetLFQ(proteinGroups_dros_hemo)

#how many missing values in protein groups dataframe
sum(is.na(proteinGroups_dros_hemo)) #28143

#remove protein groups missing more than 50% of the measurments
filtered_dros_hemo_50percent <- remove.features.50percentcuttoff(proteinGroups_dros_hemo) #[1] 6035   11
rownames(filtered_dros_hemo_50percent) <- filtered_dros_hemo_50percent$id
#keep only the complete cases
complete.case_filtered_dros_hemo <- filtered_dros_hemo_50percent[complete.cases(filtered_dros_hemo_50percent),] #[1] 4851   11


cat(paste0("Number of protein groups in Hemolymph Data before any filtering of missing measurements: ", nrow(complete.case_filtered_dros_hemo))) #4815
#proteinGroups_dros_hemo[, -which(colrows(is.na(proteinGroups_dros_hemo)) > 0.5)] #one line of code that filters by 50%
# features.missing = rowMeans(is.na(proteinGroups_dros_hemo)) 
# print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
# features.missing.50more = rownames(proteinGroups_dros_hemo)[features.missing > 0.50] 
# 
# keep.features = which(features.missing <= 0.50) 
# print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
# names(keep.features) = keep.features 
# 
# remove.features = which(features.missing > 0.50)
# print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
# names(remove.features) = remove.features


#### Subset Data with log2 values only ####

subset_hemo <- function(x){
  #subset matrix containing only abundance values and log2 transform
  hemo_log2 <- log2(as.matrix(x[,c(4:11)]))
  #set row names to match protein group identifer number
  rownames(hemo_log2) <- x$id
  colnames(hemo_log2) <- c("Ctr_1hr", "1hr", "Ctr_4hr","4hr", "Ctr_8hr", "8hr", "Ctr_24hr", "24hr")
  
  return(hemo_log2)
  }

# for imputation data
hemo_50percent_log2 <- subset_hemo(filtered_dros_hemo_50percent)


# for complete cases
hemo_complete.case_log2 <- subset_hemo(complete.case_filtered_dros_hemo)


#meta data for hemo data identifying if a sample is a control or TBI
hemo_meta <- data.frame(HemoTimePoints = colnames(hemo_50percent_log2),Sample_Type=rep(c("Control","TBI"),4))
hemo_times <- c((1/24), (1/24), (4/24), (4/24), (8/24), (8/24), (24/24), (24/24))
hemo_meta$time <- hemo_times
  
#### Histogram of data ####

#histogram of log2(intensities) for data set that has NA's
hist(hemo_50percent_log2, breaks = 20, xlab = "Log2(Hemo data with the 50% filter)", main = "Histogram of Hemo Data Before Imputation")

hist(hemo_complete.case_log2, breaks = 20, xlab = "Log2(Intensites)", main = "Histogram of Hemo Data with Complete Cases")

#### PCA prcomp function complete cases ####
pca_complete.case_hemolog2 <- prcomp(t(hemo_complete.case_log2))

# Print standard deviations of the PC's determined
print(pca_complete.case_hemolog2)

# Plot function outputs a plot with the variance of the data on the y-axis and PC on x-axis
plot(pca_complete.case_hemolog2, type = "l")

# summary function summarizes the importance of each PC
# The first line describes the standard deviation for each PC
# The second describes the proportion of the variance in the data that is explained by each component
# The third describes the cumulative proportion of the explained variance
summary(pca_complete.case_hemolog2)

#### PCA Prcomp Plot ####
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

# PCA Plot for Scores of Complete Cases
sample.colors = as.numeric(factor(hemo_meta$Sample_Type))
plot(pca_complete.case_hemolog2$x, pch = 16, col = colors[sample.colors], main = "Prcomp Scores for Hemolymph Complete Cases")
text(pca_complete.case_hemolog2$x[,1], pca_complete.case_hemolog2$x[,2], labels = colnames(hemo_complete.case_log2))
legend("topright", legend = levels(hemo_meta$Sample_Type), pch = 16, col = colors[1:length(levels(hemo_meta$Sample_Type))],
       y.intersp = 0.7)

plot(pca_complete.case_hemolog2$rotation, pch = 16, col = sample.colors, main = "Prcomp Loadings for Hemolymph Complete Cases")
text(pca_complete.case_hemolog2$rotation[,1], pca_complete.case_hemolog2$rotation[,2], labels = rownames(hemo_complete.case_log2))


#### PCA bayesian with imputation ####
sum(is.na(hemo_50percent_log2)) #2935 missing values out of 48,280 (6% of the data)

rows_missing <- data.frame(which(is.na(hemo_50percent_log2), arr.ind = TRUE))
imputed_values <- unique(rows_missing$row)
row_missing <- hemo_50percent_log2[imputed_values,]
#missing <- which(hemo_50percent_log2[is.na(hemo_50percent_log2),])

pc_hemolog2 <- pca(hemo_50percent_log2, nPcs = 3, method = "bpca") #pca method
#extract imputed data set for hemo data
imputed_hemo <- completeObs(pc_hemolog2)
missing <- imputed_hemo[is.na(hemo_50percent_log2)]

hist(missing)

# PCA Plot for Scores of Log2(hemo data)
sample.colors = as.numeric(factor(hemo_meta$Sample_Type))
plot(scores(pc_hemolog2), pch = 16, col = sample.colors, main = "PCAScores Hemolymph Log2(Protein)")
text(scores(pc_hemolog2)[,1], scores(pc_hemolog2)[,2], labels = rownames(hemo_50percent_log2), 
     col = sample.colors)
legend("topleft", legend = levels(hemo_meta$Sample_Type), pch = 16, col = 1:length(levels(hemo_meta$Sample_Type)),
       y.intersp = 0.7)

# PCA Plot the loadings of Log2(hemo data)
plot(loadings(pc_hemolog2), pch = 19, col = sample.colors, main = "PCALoadings") 
text(loadings(pc_hemolog2)[,1], loadings(pc_hemolog2)[,2], labels = colnames(hemo_50percent_log2), 
     col = sample.colors)
legend("bottomright", legend = levels(hemo_meta$Sample_Type), pch = 16, col = 1:length(levels(hemo_meta$Sample_Type)),
       y.intersp = 0.7)


#Split data up between the controls and the TBI samples
even_index <- seq(2,8,2) #TBI samples
odd_index <- seq(1,8,2) #controls

reordered_hemo_complete.case_log2 <- hemo_complete.case_log2[,c(odd_index, even_index)]

imputed_hemo_reorder <- imputed_hemo[,c(odd_index,even_index)]

fold_change.complete.cases <- reordered_hemo_complete.case_log2[,5:8] - reordered_hemo_complete.case_log2[1:4]

# subset features with fold change greater than 1.2 to cluster
important_complete.cases_hemo_foldchange <- fold_change.complete.cases[rowSums(abs(fold_change.complete.cases)>1.2)>0,]

# bin protein groups so that we describe the changes within fold changes of -2 and 2
important_complete.cases_foldchange_compress <- important_complete.cases_hemo_foldchange
important_complete.cases_foldchange_compress[important_complete.cases_foldchange_compress < -4] = -4.1
important_complete.cases_foldchange_compress[important_complete.cases_foldchange_compress > 4] = 4.1

# Protein groups with fold change greater than 1.2
table(rowSums(abs(fold_change.complete.cases)>1.2)>0)

# Protein groups with fold change greater than 2
table(rowSums(abs(fold_change.complete.cases)>2)>0)

# Protein groups with fold change greater than 4
table(rowSums(abs(fold_change.complete.cases)>4)>0)

# Protein groups with fold change greater than 8
table(rowSums(abs(fold_change.complete.cases)>8)>0)

scaleRYG <- colorRampPalette(c("red","white","#2CA8E0"), space = "rgb")(31)

pheatmap(
  mat               = reordered_hemo_complete.case_log2,
  color             = scaleRYG,
  border_color      = NA,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = c("Ctr_1hr","Ctr_4hr", "Ctr_8hr", "Ctr_24hr","1hr","4hr", "8hr", "24hr"),
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "log2 complete cases Hemolymph Data"
)


pheatmap(
  mat               = fold_change.complete.cases,
  color             = scaleRYG,
  border_color      = NA,
  #scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = c("1hr","4hr", "8hr", "24hr"),
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "log2 complete cases Hemolymph Data"
)

pheatmap(
  mat               = important_complete.cases_foldchange_compress, #[rowSums(c2$membership >0.5)>0,],
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
  main              = "Fold Change with compressed Hemo Data, binning <-4 and >4"
)
complete.case_hclust <- hclust(dist(fold_change.complete.cases, method = "euclidean"))
clusters <- cutree(complete.case_hclust, k = 20) #k = 20)

clusplot(fold_change.complete.cases, clusters, lines = 0)

plot(fold_hclust, label= FALSE)
rect.hclust(fold_hclust, k=6, border = "red")

pheatmap(
  mat               = imputed_hemo_reorder ,
  color             = scaleRYG,
  border_color      = NA,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  labels_col = c("Ctr_1hr","Ctr_4hr", "Ctr_8hr", "Ctr_24hr","1hr","4hr", "8hr", "24hr"),
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  #cellwidth = 40,
  #cellheight = .25,
  #annotation_col    = mat_col,
  #annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 8,
  main              = "log2 imputed Hemolymph Data"
)

##### pls-da complete cases #####

treatment_hemo_plsda <- plsr(as.numeric(hemo_meta$Sample_Type) ~ t(hemo_complete.case_log2), method = "oscorespls", ncomp = 4)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
sample.colors = as.numeric(factor(hemo_meta$Sample_Type))
plot(treatment_hemo_plsda$scores, pch = 19, col = colors[sample.colors])
legend("bottomright", legend = levels(factor(hemo_meta$Sample_Type)), pch = 16, col = colors[1:length(levels(factor(hemo_meta$Sample_Type)))])
text(treatment_hemo_plsda$scores[,1], treatment_hemo_plsda$scores[,2], as.character(hemo_meta$HemoTimePoints), pos = 4, offset = 0.5, cex = 0.8)

par(mar = c(4,6,4,1), las  = 1, mgp = c(2.5,0.5,0), tcl =  -0.3, ps = 12)
dotchart(rev(treatment_hemo_plsda$loadings[order(-abs(treatment_hemo_plsda$loadings[,1])),1][1:25]),
         pch=19,
         color = "#5C66AF",
         xlab = "Component 1 Loadings")



###### time series #######

#randomly sample from the hemo_complete.case file
set.seed(123)
sample.series <- data.frame(hemo_complete.case[sample(nrow(hemo_complete.case),100),])
proteins.id <- rownames(sample.series)
sample.series$protein.ID <- proteins.id

#Transform data into wide form, using protein.ID as the id.variable
sample.series.wide <- melt(sample.series, 'protein.ID')
colnames(sample.series.wide) <- c("protein.ID", "sample.ID", "Intensity")
sample.series.wide$Log2.Intensity <- log2(sample.series.wide$Intensity)
sample.series.wide$time <- c(rep(1,200), rep(4, 200), rep(8,200), rep(24, 200))
sample.series.wide$type <- rep(c(rep("Control", 100), rep("TBI", 100)),4)


mod_sample.ID <- lm(Log2.Intensity~sample.ID, sample.series.wide)
hist(mod_sample.ID$residuals)
mean(mod_sample.ID$residuals)

mod_time <- lm(Log2.Intensity~time, sample.series.wide)
hist(mod_time$residuals)

mod_sample.ID_time <- lm(Log2.Intensity~sample.ID+time, sample.series.wide)
hist(mod_sample.ID_time$residuals)

mod_time_sample.ID <- lm(Log2.Intensity~time+sample.ID, sample.series.wide)
hist(mod_time_sample.ID$residuals)


par(mfrow = c(2,2)) # set 2 rows and 2 column plot layout
plot(mod_sample.ID)


ranked_pls_loading_hemo <- as.numeric(sub(".*log2)", "", names(rev(treatment_hemo_plsda$loadings[order(-abs(treatment_hemo_plsda$loadings[,1])),1][1:50]))))
labels_hemo_hours <- c("1","4","8","24")

hemo_timeSeries_protein_plot <- function(n){
  even_timepoints <- hemo_complete.case_log2[which(rownames(hemo_complete.case_log2)==n),seq(2,8,2)]
  odd_timepoints <- hemo_complete.case_log2[which(rownames(hemo_complete.case_log2)==n),seq(1,7,2)]
  
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
  plot(hemo_times[seq(2,8,2)], even_timepoints,
       ylim = c(minimum_y_Intensity,max_y_Intensity),
       xlab = "Time (Hours)", ylab = paste("Log2(Intensity of ",n,")"),  
       cex=2, 
       type = "b,c",
       pch = 19,
       lty = 2,
       main = paste(n),
       xaxt = "n")
  lines(hemo_times[seq(2,8,2)],odd_timepoints, col = "red") 
  axis(1,hemo_times[seq(2,8,2)], labels = labels_hemo_hours)
  legend("topright", legend = levels(factor(head_meta$Sample_Type)), 
         lty = 1, col = c("red","black"))
}



n = 70
hemo_timeSeries_protein_plot(n)
proteins_hemo <- proteinGroups_dros_hemo[which(proteinGroups_dros_hemo$id==n),]

proteins_hemo_n <- rownames(proteinGroups_dros_hemo[grepl("Q8IN44", proteinGroups_dros_hemo$Protein.IDs),])
df_proteins_hemo <- proteinGroups_dros_hemo[proteins_hemo_n,]
n <- as.numeric(df_proteins_hemo$id)

top_proteins_hemo_loadings <- rev(treatment_hemo_plsda$loadings[order(-abs(treatment_hemo_plsda$loadings[,1])),1][1:50])
names(top_proteins_hemo_loadings) <- ranked_pls_loading_hemo

df_hemo <- data.frame(id = as.numeric(), Anova.p.Samply_Type = as.numeric(), 
                      Anova.p.time = as.numeric(), Anova.p.Samply_Type.Time = as.numeric(), na = as.character(), Loading = as.numeric())
for(i in 1:length(ranked_pls_loading_hemo)){
  n <- ranked_pls_loading_hemo[i]
  
  aov_n_hemo <- aov(hemo_complete.case_log2[which(rownames(hemo_complete.case_log2)==n),] ~ hemo_meta$Sample_Type*hemo_meta$time)
  
  df_hemo[i,] <- c(n,summary(aov_n_hemo)[[1]][["Pr(>F)"]],top_proteins_hemo_loadings[as.character(n)][[1]])
}

matching_id_hemo <- which(complete.case_filtered_dros_hemo$id %in% df_hemo$id)
matching_proteins_genes_hemo <- complete.case_filtered_dros_hemo[matching_id_hemo, c(1:3)]

pls_top_proteins_hemo <- merge(matching_proteins_genes_hemo, df_hemo, by.y = "id")

rownames(pls_top_proteins_hemo) <- pls_top_proteins_hemo$id
pls_top_proteins_hemo <- merge(pls_top_proteins_hemo, fold_change.complete.cases, by = "row.names")

write.csv(pls_top_proteins_hemo, file = "G:/Projects/Proteomics/DorsophilaHead_Experiment/Routput/Hemo/20190325pls_top_protein_gene_id_aov_HEMO_completeCase.csv", row.names = FALSE)


for(i in 1:length(ranked_pls_loading_hemo)){
  n = pls_top_proteins_hemo[i,1]
  
  pdf(paste("G:/Projects/Proteomics/DorsophilaHead_Experiment/Figures/Hemo/",n,"_timeseriesplot.pdf"), useDingbats = FALSE)
  hemo_timeSeries_protein_plot(n)
  dev.off()
}
