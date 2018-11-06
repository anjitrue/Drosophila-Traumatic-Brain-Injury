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

Ensembl_protein_id <- proteinGroups_Ensembl_reducedtoLFQ$`Protein IDs`
Uniprot_protein_id <- proteinGroups_Uniprot_reducedtoLFQ$`Protein IDs`


Ensembl_join <- data.frame(proteinGroups_Ensembl_reducedtoLFQ[,c(1,2)])
rownames(Ensembl_join) <- Ensembl_join$`Protein.IDs`
#protein_id <- Uniprot_test$Protein.IDs

#--------------------------------------------- BiomaRt Code -------------------------------------------
ensembl = useMart("ensembl") # the useMart() can connect to a specificed MioMart database
Ensembl_datasets <-listDatasets(ensembl) #view datasets available in Ensembl

ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl) #select dataset mouse Mouse genes (GRCm38.p6)

filters = listFilters(ensembl) #filter options
attributes <- listAttributes(ensembl) #values interested in retreiving from query

searchFilters(mart = ensembl, pattern = "mmus")
searchAttributes(mart = ensembl, pattern = "mmus")

#filterType("with_uniprotswissprot", ensembl)

transcript_mouse_ID_Ensembl <- getBM(attributes = c("ensembl_transcript_id", "uniprotsptrembl", "uniprotswissprot"),
                             values = Ensemble_protein_id,
                             #filters = "gene_name",
                             mart = ensembl)

transcript_mouse_ID_Ensembl <- data.frame(transcript_mouse_ID_Ensembl)
colnames(transcript_mouse_ID_Ensembl)[which(names(transcript_mouse_ID_Ensembl) == "ensembl_transcript_id")] <- "Protein.IDs"


Ensembl_join = right_join(transcript_mouse_ID_Ensembl, Ensembl_join, by = "Protein.IDs")

write.table(transcript_mouse_ID_Ensembl, "E:/Projects/Proteomics/DiversityOutcross/txt_proteinEnsembl/Rmodified_EnsembleID_matchingTable.txt", sep="\t")


#--------------------------------------------- AnnotationHub Code -------------------------------------------                             

hub <- query(AnnotationHub(), c('ensembl', 'gtf', 'mus musculus'))
hub <- hub[grep('GRCm38', hub$title)]
anno_hub <- hub[["AH60127"]]

#union <- Ensembl_test[which(rownames(filtered) %in% Ensembl_test$Protein.IDs),]