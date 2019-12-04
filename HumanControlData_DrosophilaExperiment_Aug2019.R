##### Load data #####
HC_drosophila <- read.csv("C:/Users/etrujillo2/Documents/Projects/Proteomics/FDR summary_Drosophila_HC_20190807.csv",header = TRUE, sep = ",", stringsAsFactors = FALSE)

HC_drosophila[,1] <- c(1, 2, 3, 4, 5)

m <- ggplot(HC_drosophila, aes(x = CSV.File, y = Total.MS.MS.Spectra, group = 1)) +
  #scale_x_continuous(breaks = seq(100, 2000, 200)) +
  ylim(235000, 247000)+
  scale_y_continuous(breaks = seq(23500,247000, 1000)) +
  geom_point(size = 5, color = "mediumorchid2") +
  geom_line(linetype = "dashed")+
  theme_light()
m + labs(title = "Number of MS/MS", x = "Date", y = "Number of MS/MS")

p <- ggplot(HC_drosophila, aes(x = CSV.File, y = X.1, group = 1)) +
  geom_bar(stat = "identity", color = "#6C686E" , fill = "#89469B", position = position_dodge())+
  theme_light()
p + labs(title = "Percent of MS/MS Reduced", x = "Date", y = "Percent MS/MS")

m <- ggplot(HC_drosophila, aes(x = CSV.File, y = PSMs, group = 1)) +
  #scale_x_continuous(breaks = seq(100, 2000, 200)) +
  #ylim(235000, 247000)+
  scale_y_continuous(breaks = seq(92000,102400, 2000)) +
  geom_point(size = 5, color = "mediumorchid2") +
  #geom_point(data = reformated_rapid_numPeptides_2000, size = 5, color = "olivedrab3")+
  #geom_point(data = reformated_normal_numPeptides_2000, size = 5, color = "steelblue2")+
  geom_line(linetype = "dashed")+
  #geom_line(data = reformated_rapid_numPeptides_2000, linetype = "dashed")+
  #geom_line(data = reformated_normal_numPeptides_2000, linetype = "dashed")+
  theme_light()
m + labs(title = "Number of PSMs", x = "Date", y = "Number of PSMs")

p <- ggplot(HC_drosophila, aes(x = CSV.File, y = X.4, group = 1)) +
  geom_bar(stat = "identity", color = "#6C686E" , fill = "#89469B", position = position_dodge())+
  theme_light()
p + labs(title = "Percent of PSMs Reduced", x = "Date", y = "Percent PSMs")

m <- ggplot(HC_drosophila, aes(x = CSV.File, y = Peptides, group = 1)) +
  #scale_x_continuous(breaks = seq(100, 2000, 200)) +
  #ylim(235000, 247000)+
  scale_y_continuous(breaks = seq(62000,67300, 500)) +
  geom_point(size = 5, color = "mediumorchid2") +
  #geom_point(data = reformated_rapid_numPeptides_2000, size = 5, color = "olivedrab3")+
  #geom_point(data = reformated_normal_numPeptides_2000, size = 5, color = "steelblue2")+
  geom_line(linetype = "dashed")+
  #geom_line(data = reformated_rapid_numPeptides_2000, linetype = "dashed")+
  #geom_line(data = reformated_normal_numPeptides_2000, linetype = "dashed")+
  theme_light()
m + labs(title = "Number of Peptides", x = "Date", y = "Number of Peptides")

p <- ggplot(HC_drosophila, aes(x = CSV.File, y = X.7, group = 1)) +
  geom_bar(stat = "identity", color = "#6C686E" , fill = "#89469B", position = position_dodge())+
  theme_light()
p + labs(title = "Percent of Peptides Reduced", x = "Date", y = "Percent Peptides")
