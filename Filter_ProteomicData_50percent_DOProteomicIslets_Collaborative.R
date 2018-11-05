#--------------------------------------------- Functions ---------------------------------------------

#create a function to remove protein groups that contain less than 50% of data
remove.features <- function (x) {
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

#--------------------------------------------- Installation of Packages -------------------------------
library(readxl)

# Upload data from absolute path
proteinGroups_Ensembl_reducedtoLFQ <- read_excel("E:/Projects/Proteomics/DiversityOutcross/txt_proteinEnsembl/proteinGroups_Ensmbl_reducedtoLFAQ.xlsx") 
Ensembl_test <- data.frame(proteinGroups_Ensembl_reducedtoLFQ[,-c(2:18,31:32)]) #subset columns that contain LFQ plus protein group column, will have to change columns depending on dataset

rownames(Ensembl_test) <- Ensembl_test$Protein.IDs #name row names by protein group

#--------------------------------------------- Data Manipulation ---------------------------------------

# Convert 0 as NA and sum the number of NA present in data base
Ensembl_test[Ensembl_test == 0] <- NA 
sum(is.na(Ensembl_test)) #66285 uniprot 56254 Ensembl

# Filter data to remove protein groups that contain less than 50% of measurements
filtered <- remove.features(Ensembl_test)

#--------------------------------------------- Write to txt file ---------------------------------------

write.table(filtered, "E:/Projects/Proteomics/DiversityOutcross/txt_proteinEnsembl/Rmodified_50percentfilterd_Ensembl_test_DOMiceProteomics.txt", sep="\t")




