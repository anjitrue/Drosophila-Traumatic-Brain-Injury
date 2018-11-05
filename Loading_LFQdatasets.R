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


#create a vector with all the time points 
timePoint <- c("Immediate following HIT","30 min","1 hour","2 hour","4 hour","6 hour","8 hour","12 hour","16 hour","24 hour",
               "48 hour","72 hour","7 days","14 days","21 days","28 days","35 days")

#Heads_1_hour_Branch_Heads_Proteomics_Fractions_ <- read_delim("E:/Projects/Proteomics/DorsophilaHead_Experiment/LFQ/2Branches/Heads - 1 hour (Branch_ Heads Proteomics + Fractions).txt", 
                                                #"\t", escape_double = FALSE, trim_ws = TRUE)
#set working directory
getwd()
setwd("E:/Projects/Proteomics/DorsophilaHead_Experiment/LFQ/2Branches")

# create a vecotr of files in directory with .txt ending, read files in with read.delim
temp = list.files(pattern = "*.txt")
myfiles = lapply(temp, read.delim)

#populate a list with the correct file name and corresponding database (create an lapply function)
filename = list()
for(i in 1:length(temp)){ 
  filename[i] <- temp[i]
}

names(filename) <- filename

for(i in 1:length(filename)){ 
  filename[i] <- myfiles[i]
}

filename_totaldf <- data.frame()


y <- merge(filename[[1]], filename[[2]], by = 1, all = TRUE)


###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

setwd("E:/Projects/Proteomics/DorsophilaHead_Experiment/")


#load data
proteinGroups <- read_delim("proteinGroups_EAT_DrosHeadsHemo_201801025.txt","\t", escape_double = FALSE, trim_ws = TRUE)
sapply(strsplit(proteinGroups, ";"), "[", 1)                                                      

rownames(proteinGroups) <- proteinGroups$id

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


remove.features <- function (x)
  {
    x = log2protein_heads
    features.missing = rowMeans(is.na(x)) 
    print(sum(features.missing > 0.50)) # 1524 protein groups
    features.missing.50more = rownames(x)[features.missing > 0.50] 
    
    keep.features = which(features.missing < 0.50) #7568
    print(paste0("Remaining Features: ", length(keep.features)))
    keep.features = names(keep.features) 
    
    remove.features = which(features.missing > 0.50) #1524
    print(paste0("Features Removed: ", length(remove.features)))
    remove.features = names(remove.features)
    
    filtered = x[-which(rownames(x) %in% remove.features),]
    return(filtered)
}
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

