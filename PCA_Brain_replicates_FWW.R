colors <- c("#206F94", # teal
            "#F47F72", #coral
            "#75C69D", #baby green
            "#5C66AF", #light purpleblue
            "#2CA8E0", #coon blue
            "#1F6F94", #darker blue
            "#ED237A", #hot pink)
            "#5C66AF", #light purpleblue
            "#2A4093", #dark purpleblue
            "#2C9379", #dark baby green
            "#83D5F7", #light coon blue
            "#93211E", #dark red
            "#E73C25", #red
            "#81143A" #dark pink
            ) 

#PCA of all the samples
sample.colors = as.numeric(factor(brain_FWW_meta$Condition))
plot(pca_complete_LFQ_FWW$x, pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Brain TBI on Food & Water or Water", 
     ylim = c(-20, 20), xlim = c(-25, 15), 
     xlab = "PC1 17.95%", ylab ="PC2 9.78%")
#text(pca_complete_LFQ_FWW$x[,1], pca_complete_LFQ_FWW$x[,2], labels = brain_FWW_meta$Samples)
legend("topleft", legend = levels(factor(brain_FWW_meta$Condition)), pch = 16, 
       col = colors[1:length(levels(factor(brain_FWW_meta$Condition)))], y.intersp = 0.7)

#
sample.colors = as.numeric(factor(brain_FWW_meta[58:97,]$replicates))
plot(pca_complete_LFQ_FWW$x[58:97,1:2], pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Brain TBI  Water", 
     ylim = c(-20, 20), xlim = c(-20, -10), 
     xlab = "PC1\n32.48%", ylab ="PC2\n10.27%")
text(pca_complete_LFQ_FWW$x[,1][58:97], pca_complete_LFQ_FWW$x[,2][58:97], 
     labels = brain_FWW_meta$Samples[58:97])
legend("topright", legend = levels(factor(brain_FWW_meta$replicates)), pch = 16, 
       col = colors[1:length(levels(factor(brain_FWW_meta$replicates)))], y.intersp = 0.7)


sample.colors = as.numeric(factor(brain_FWW_meta$Sample_Type[7:12]))
plot(pca_FW$x[7:12,1:2], pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Complete Head Data Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 20), 
     xlab = "PC1\n32.48%", ylab ="PC2\n10.27%")
text(pca_FW$x[,1][7:12], pca_FW$x[,2][7:12], 
     labels = colnames(LFQ_FWW_Ordered_LOG2)[7:12])
legend("topright", legend = levels(factor(brain_FWW_meta$Sample_Type[1:57])), pch = 16, 
       col = colors[1:length(levels(factor(brain_FWW_meta$Sample_Type[1:57])))], y.intersp = 0.7)

sample.colors = as.numeric(factor(brain_FWW_meta$Sample_Type[7:12]))
plot(pca_FW$x[7:12,1:2], pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Complete Head Data Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 20), 
     xlab = "PC1\n32.48%", ylab ="PC2\n10.27%")
text(pca_FW$x[,1][7:12], pca_FW$x[,2][7:12], 
     labels = colnames(LFQ_FWW_Ordered_LOG2)[7:12])
legend("topright", legend = levels(factor(brain_FWW_meta$Sample_Type[1:57])), pch = 16, 
       col = colors[1:length(levels(factor(brain_FWW_meta$Sample_Type[1:57])))], y.intersp = 0.7)

# PCA by replicate, shows batch affect by replicate
sample.colors = as.numeric(factor(brain_FWW_meta$replicates[1:57]))
plot(pca_FW$x[1:57,1:2], pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 20), 
     xlab = "PC1 24.17%", ylab ="PC2 10.23%")
#text(pca_FW$x[,1][1:57], pca_FW$x[,2][1:57], labels = colnames(LFQ_FWW_Ordered_LOG2)[1:57])
legend("topright", legend = levels(factor(brain_FWW_meta$replicates[1:57])), pch = 16, 
       col = colors[3:5][1:length(levels(factor(brain_FWW_meta$replicates[1:57])))], y.intersp = 0.7)


#Sample_Type TBI vs Control FW
sample.colors = as.numeric(factor(brain_FWW_meta$hour[1:57]))
plot(pca_FW$x[1:57,1:2], pch = 19, cex = 2,
     col = colors[sample.colors], 
     main = "Principle Component Analysis\n Brain TBI on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 20), 
     xlab = "PC1 24.17%", ylab ="PC2 10.23%")
#text(pca_FW$x[,1][1:57], pca_FW$x[,2][1:57], labels = colnames(LFQ_FWW_Ordered_LOG2)[1:57])
legend("topleft", legend = levels(factor(brain_FWW_meta$hour[1:57])), pch = 16, 
       col = colors[1:length(levels(factor(brain_FWW_meta$hour[1:57])))], y.intersp = 0.7)


#replicate FW
sample.colors = as.numeric(factor(brain_FWW_meta$replicates[1:57]))
plot(pca_FW$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 20), 
     xlab = "PC1 24.17%", ylab ="PC2 10.23%")
#text(pca_FW$x[,1][1:57], pca_FW$x[,2][1:57], labels = colnames(LFQ_FWW_Ordered_LOG2)[1:57])
legend("topright", legend = levels(factor(brain_FWW_meta$replicates[1:57])), pch = 16, 
       col = colors[3:5][1:length(levels(factor(brain_FWW_meta$replicates[1:57])))], y.intersp = 0.7)



sample.colors = as.numeric(factor(brain_FWW_meta$replicates[58:97]))
plot(pca_W$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI on Water Only", 
     ylim = c(-20, 20), xlim = c(-25, 20), 
     xlab = "PC1 16.43%", ylab ="PC2 11.66%")
text(pca_W$x[,1], pca_W$x[,2], 
     labels = colnames(LFQ_FWW_Ordered_LOG2)[58:97])
legend("topright", legend = levels(factor(brain_FWW_meta$replicates[58:97])), pch = 16, 
       col = colors[3:5][1:length(levels(factor(brain_FWW_meta$replicates[58:97])))], y.intersp = 0.7)


sample.colors = as.numeric(factor(mean_normaized_FW_byProtein_sepTBIandCTRL$replicates))
plot(pca_FW_meanNormalized$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Normalized to Mean Protein Intensity\n on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 16.18%", ylab ="PC2 13.99%")
#text(pca_FW$x[,1][1:57], pca_FW$x[,2][1:57], labels = colnames(LFQ_FWW_Ordered_LOG2)[1:57])
legend("topright", legend = levels(factor(mean_normaized_FW_byProtein_sepTBIandCTRL$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(mean_normaized_FW_byProtein_sepTBIandCTRL$replicates)))], y.intersp = 0.7)



# PCA of the samples normalized to combine TBI and Ctrl time 0. PCA is looking at Odd samples or control samples
sample.colors = as.numeric(factor(new_FW_Log2_odd$replicates))
plot(pca_FW_meanNormalized_odd$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Normalized to Mean Protein Intensity of all time 0\n Control Samples on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 25.23%", ylab ="PC2 20.76%")
text(pca_FW_meanNormalized_odd$x[,1], pca_FW_meanNormalized_odd$x[,2], labels = new_FW_Log2_odd$Samples)
legend("topright", legend = levels(factor(new_FW_Log2_odd$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(new_FW_Log2_odd$replicates)))], y.intersp = 0.7)

sample.colors = as.numeric(factor(new_FW_Log2_even$replicates))
plot(pca_FW_meanNormalized_even$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Normalized to Mean Protein Intensity\n TBI Samples on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 27.87%", ylab ="PC2 16.24%")
#text(pca_FW$x[,1][1:57], pca_FW$x[,2][1:57], labels = colnames(LFQ_FWW_Ordered_LOG2)[1:57])
legend("topright", legend = levels(factor(new_FW_Log2_even$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(new_FW_Log2_even$replicates)))], y.intersp = 0.7)

sample.colors = as.numeric(factor(new_FW_Log2$replicates))
plot(pca_FW_meanNormalized_combined$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Normalized to Mean Protein Intensity of all time 0\n All Samples on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 16.13%", ylab ="PC2 13.18%")
text(pca_FW_meanNormalized_combined$x[,1], pca_FW_meanNormalized_combined$x[,2], labels = new_FW_Log2$Samples)
legend("topright", legend = levels(factor(new_FW_Log2$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(new_FW_Log2$replicates)))], y.intersp = 0.7)


ggplot(mean_normaized_FW_byProtein_sepTBIandCTRL, aes(x= mean_normaized_FW_byProtein_sepTBIandCTRL$`4`)) + 
        geom_histogram(binwidth = .05)
ggplot(new_FW_Log2, aes(x=new_FW_Log2$`4`)) + geom_histogram(binwidth = .05)
ggplot(t_LFQ_FW_ordered_LOG2, aes(x=t_LFQ_FW_ordered_LOG2$`4`)) + geom_histogram(binwidth = .05)

variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)


sample.colors = as.numeric(factor(new_FW_odd$replicates))
plot(pca_FW_TIC_taverage$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Normalized to TIC time averaged \n Control Samples on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 16.13%", ylab ="PC2 13.18%")
text(pca_FW_TIC_taverage$x[,1], pca_FW_TIC_taverage$x[,2], labels = new_FW_odd$Samples)
legend("topright", legend = levels(factor(new_FW_odd$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(new_FW_odd$replicates)))], y.intersp = 0.7)

sample.colors = as.numeric(factor(t_LFQ_FW_ordered_lessThan100CV$replicates))
plot(pca_FW_lessThan100CV$x, pch = 19, cex = 2,
     col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Less than 100% CV\n Control Samples on Food & Water", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 16.13%", ylab ="PC2 13.18%")
text(pca_FW_lessThan100CV$x[,1], pca_FW_lessThan100CV$x[,2], labels = t_LFQ_FW_ordered_lessThan100CV$Samples)
legend("topright", legend = levels(factor(t_LFQ_FW_ordered_lessThan100CV$replicates)), pch = 16, 
       col = colors[3:5][1:length(levels(factor(t_LFQ_FW_ordered_lessThan100CV$replicates)))], y.intersp = 0.7)

# Loadings plot of samples with proteins less than 100% CV
plot(pca_FW_lessThan100CV$rotation, pch = 19, cex = 2,
     #col = colors[3:5][sample.colors], 
     main = "Principle Component Analysis\n Brain TBI Less than 100% CV\n Control Samples on Food & Water", 
     #ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 25.43%", ylab ="PC2 9.13%")
text(pca_FW_lessThan100CV$rotation[,1][which(pca_FW_lessThan100CV$rotation[,1] < -0.045)], pca_FW_lessThan100CV$rotation[,2][which(pca_FW_lessThan100CV$rotation[,1] < -0.045)], 
     col = "red", labels = names(pca_FW_lessThan100CV$rotation[,1][which(pca_FW_lessThan100CV$rotation[,1] < -0.045)]))
text(pca_FW_lessThan100CV$rotation[,1][which(pca_FW_lessThan100CV$rotation[,1] > 0.045)], pca_FW_lessThan100CV$rotation[,2][which(pca_FW_lessThan100CV$rotation[,1] > 0.045)], 
     col = "red", labels = names(pca_FW_lessThan100CV$rotation[,1][which(pca_FW_lessThan100CV$rotation[,1] > 0.045)]))

plot(pca_complete_HC$x, pch = 19, cex = 2,
     col = colors[1:7], 
     main = "Principle Component Analysis\n Human Controls LFQ", 
     ylim = c(-20, 20), xlim = c(-20, 25), 
     xlab = "PC1 36.26%", ylab ="PC2 26.63%")
text(pca_complete_HC$x[,1], pca_complete_HC$x[,2], labels = colnames(hc_log2))







ggplot(t_LFQ_FW_ordered_LOG2_even , aes(x = hour , 
                                                      y = t_LFQ_FW_ordered_LOG2_even$`4`, fill = replicates)) +
        geom_boxplot()

ggplot(mean_normaized_FW_byProtein_sepTBIandCTRL[1:55,], aes(x = hour , 
                                                      y = mean_normaized_FW_byProtein_sepTBIandCTRL[1:55,]$`4`, fill = replicates)) +
        geom_boxplot()


ggplot(new_FW_Log2_even, aes(x = hour , 
                                                      y = new_FW_Log2_even$`4`, fill = replicates)) +
        geom_boxplot()



ggplot(t_LFQ_FW_ordered , aes(x = Samples , 
                             y = log2(TIC), fill = hour)) +
        geom_boxplot()

plot(log2(new_FW_odd$TIC))
text(log2(new_FW_odd$TIC), labels = t_LFQ_FW_ordered_odd$Samples)


plot(log2(col_sums), main = "TIC - Time point", ylab = "TIC", xlab = "Samples")
text(log2(col_sums), labels = colnames(LFQ_FWW_Ordered[,1:57]))

plot(log2(new_FW_odd$TIC), main = "TIC - Time point", ylab = "TIC", xlab = "Samples")
text(log2(new_FW_odd$TIC), labels = new_FW_odd$Samples)


hist(cv_time.point.matrix)
hist(hc$CV)

h <- ggplot(cv_time.point.data.frame, aes(x = `1`)) + geom_histogram(binwidth = 1, color = "black", fill = "black", alpha = 0.2)
# h + geom_histogram(aes(y=..density..), position="identity", alpha=0.5) # not working

h + geom_histogram(aes(x = `2`), binwidth = 1, colour="black", fill="#206F94", alpha = 0.2) + 
        geom_histogram(aes(x = `3`), binwidth = 1, colour="black", fill="#F47F72", alpha = 0.2)

d <- ggplot(cv_time.point.data.frame, aes(x = `1`)) + geom_density(fill = "black", alpha = 0.2)
d + geom_density(aes(x = `3`), fill = "#206F94", alpha = 0.2) + geom_density(aes(x = `5`), fill = "#F47F72", alpha = 0.2) + 
        geom_density(aes(x = `7`), fill = "#75C69D", alpha = 0.2) + geom_density(aes(x = `9`), fill = "#5C66AF", alpha = 0.2) +
        geom_density(aes(x = `11`), fill = "#ED237A", alpha = 0.2) + geom_density(aes(x = `13`), fill = "#5C66AF", alpha = 0.2) +
        geom_density(aes(x = `15`), fill = "#2A4093", alpha = 0.2) + geom_density(aes(x = `17`), fill = "#83D5F7", alpha = 0.2) +
        geom_density(aes(x = `19`), fill = "#2C9379", alpha = 0.2) 
        
plot(cv_time.point.data.frame$`1`, cv_time.point.data.frame$`3`)
text(cv_time.point.data.frame$`1`, cv_time.point.data.frame$`3`, 
     labels = rownames(cv_time.point.data.frame[which(cv_time.point.data.frame[,c(1,3)] > 100),c(1,3)]))

rownames(cv_time.point.data.frame[which(cv_time.point.data.frame[,c(1,3)] > 100),])

cv_100_meta <- proteinGroups_dros_Brain_FWW_extended[which(proteinGroups_dros_Brain_FWW_extended$id %in%cv_100_protein.id),]

plot(table(cv_100$sample))
