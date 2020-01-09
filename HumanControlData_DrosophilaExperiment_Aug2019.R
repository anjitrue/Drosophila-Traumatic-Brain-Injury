##Packages
library(ggplot2)

##### Load data #####
HC_drosophila <- read.csv("C:/Users/etrujillo2/Documents/Projects/Proteomics/FDR summary_Drosophila_HC_20190807.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)

#rename the first column using numeric sequence 1:5
HC_drosophila[,1] <- c(1, 2, 3, 4, 5)

HC_drosophila$MS2Percentage <- (HC_drosophila$Total.MS.MS.Spectra/HC_drosophila$Total.MS.MS.Spectra[1])*100
HC_drosophila$MS2overPSMs <- (HC_drosophila$PSMs/HC_drosophila$Total.MS.MS.Spectra)*100
HC_drosophila$PeptidesoverPSMS <- (HC_drosophila$Peptides/HC_drosophila$PSMs)*100

# Plot the number of MS2 scans relative the first day of analysis
m <- ggplot(HC_drosophila, aes(x = CSV.File, y = MS2Percentage, group = 1)) +
  scale_y_continuous(breaks = seq(0,100, 10), limits = c(0, 100)) +
  geom_point(size = 5, color = "mediumorchid2") +
  geom_line(linetype = "dashed")+
  theme_light()
m + labs(title = "Normalized MS/MS scans relative to first day of analysis", x = "Human Contorl", y = "Percent of MS/MS Scans")

# Another way of plotting the loss in performance is to plot the loss of MS/MS scans relative to the first day of analysis
p <- ggplot(HC_drosophila, aes(x = CSV.File, y = X.1, group = 1)) +
  scale_y_continuous(breaks = seq(0,10, 1), limits = c(0, 10)) +
  geom_point(size = 5, color = "mediumorchid2") +
  geom_line(linetype = "dashed")+
  theme_light()
p + labs(title = "Percent of MS/MS scans reduced over 2 week analysis time", x = "Human Control", y = "Percent of MS/MS Scans")

# Look at percent of PSMs relative to the number of MS/MS scans "identification rate"
m <- ggplot(HC_drosophila, aes(x = CSV.File, y = MS2overPSMs, group = 1)) +
  scale_y_continuous(breaks = seq(0,100, 10), limits = c(0, 100)) +
  geom_bar(stat = "identity", color = "mediumorchid2")+
  theme_light()
m + labs(title = "Percent of PSMs relative to the number of MS/MS scans", x = "Human Control", y = "Percent")

# Calculate the percent of unique peptides relative to PSMs "unique score"
m <- ggplot(HC_drosophila, aes(x = CSV.File, y = PeptidesoverPSMS, group = 1)) +
  scale_y_continuous(breaks = seq(0,100, 10), limits = c(0, 100)) +
  geom_bar(stat = "identity", color = "mediumorchid2")+
  theme_light()
m + labs(title = "Number of unique peptides relative to number of PSMs", x = "Human Control", y = "Percentage of unique peptides relative to PSMs")

# p <- ggplot(HC_drosophila, aes(x = CSV.File, y = X.7, group = 1)) +
#   geom_bar(stat = "identity", color = "#6C686E" , fill = "#89469B", position = position_dodge())+
#   theme_light()
# p + labs(title = "Percent of Peptides Reduced", x = "Date", y = "Percent Peptides")
