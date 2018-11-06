#--------------------------------------------- Installation of Packages -------------------------------
install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
BiocManager::install("AnnotationHub", version = "3.8")
BiocManager::install("rtracklayer", version = "3.8")

library(readxl)
library(BiocManager)
library(biomaRt)
library(AnnotationHub)
library(rtracklayer)
library(dplyr)

#--------------------------------------------- Data Upload --------------------------------------------
proteinGroups_Ensembl_reducedtoLFQ <- read_excel("E:/Projects/Proteomics/DiversityOutcross/txt_proteinEnsembl/proteinGroups_Ensmbl_reducedtoLFAQ.xlsx") 
proteinGroups_Uniprot_reducedtoLFQ <- read_excel("E:/Projects/Proteomics/DiversityOutcross/txt_Uniprot/proteinGroups_Uniprot_reducedtoLFQ.xlsx") 

# vectors containing the protein.ids
Ensembl_protein_id <- proteinGroups_Ensembl_reducedtoLFQ$`Protein IDs`
Uniprot_protein_id <- proteinGroups_Uniprot_reducedtoLFQ$`Protein IDs`

# Subset dataframe to only the protein IDs
Uniprot_join <- data.frame(proteinGroups_Uniprot_reducedtoLFQ[,c(1,2)])
rownames(Uniprot_join) <- Uniprot_join$`Protein.IDs`
#protein_id <- Uniprot_test$Protein.IDs

#--------------------------------------------- BiomaRt Code -------------------------------------------
ensembl = useMart("ensembl") # the useMart() can connect to a specificed MioMart database
Ensembl_datasets <-listDatasets(ensembl) #view datasets available in Ensembl

ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl) #select dataset mouse Mouse genes (GRCm38.p6)

filters = listFilters(ensembl) #filter options
attributes <- listAttributes(ensembl) #values interested in retreiving from query

searchFilters(mart = ensembl, pattern = "uniprot")
searchAttributes(mart = ensembl, pattern = "ensembl")

# Create matching table for Uniprot protein.ids
transcript_mouse_ID_Ensembl <- getBM(attributes = c("ensembl_transcript_id","uniprotsptrembl", "uniprotswissprot"),
                             values = Ensemble_protein_id,
                             filters = "ensembl_transcript_id",
                             mart = ensembl)

transcript_mouse_ID_Uniprot <- getBM(attributes = c("ensembl_transcript_id", "uniprotsptrembl", "uniprotswissprot"),
                                     values = Uniprot_protein_id,
                                     #filters = "uniprotswissprot", # important filter to obtain only UniprotIDs in list
                                     mart = ensembl)

# Matching to existing list
transcript_mouse_ID_Ensembl <- data.frame(transcript_mouse_ID_Ensembl)
colnames(transcript_mouse_ID_Ensembl)[which(names(transcript_mouse_ID_Ensembl) == "ensembl_transcript_id")] <- "Protein.IDs"


Uniprot_join = right_join(transcript_mouse_ID_Ensembl, Uniprot_join, by = "Protein.IDs")

write.table(transcript_mouse_ID_Uniprot, "E:/Projects/Proteomics/DiversityOutcross/txt_Uniprot/Rmodified_Uniprot_matchingTable.txt", sep="\t")

i = 1
match_id <- NULL
for(i in 1:nrow(transcript_mouse_ID_Uniprot)){
  if(!(transcript_mouse_ID_Uniprot[i,2] == ""| is.na(transcript_mouse_ID_Uniprot[i,2]) )){
   if (length(grep(transcript_mouse_ID_Uniprot[i,2], Uniprot_protein_id) > 0)){ 
     x = grep(transcript_mouse_ID_Uniprot[i,2], Uniprot_protein_id)
     } else(x = NA)
  } else(x = NA)
  if(!(transcript_mouse_ID_Uniprot[i,3] == "" | is.na(transcript_mouse_ID_Uniprot[i,2]))){
    if (length(grep(transcript_mouse_ID_Uniprot[i,3], Uniprot_protein_id)) > 0){ 
      y = grep(transcript_mouse_ID_Uniprot[i,3], Uniprot_protein_id) } else(y = NA)
  } else(y = NA)
  match_id <- append(match_id, paste(x,y))
}
match_id <- gsub(" ", "", match_id)
match_id <- gsub("NA", "", match_id)
match_id <- sub(";", "", match_id)

table(match_id)

#--------------------------------------------- AnnotationHub Code -------------------------------------------                             

hub <- query(AnnotationHub(), c('ensembl', 'gtf', 'mus musculus'))
hub <- hub[grep('GRCm38', hub$title)]
anno_hub <- hub[["AH60127"]]

#--------------------------------------------- Data Upload --------------------------------------------

#union <- Ensembl_test[which(rownames(filtered) %in% Ensembl_test$Protein.IDs),]