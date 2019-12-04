library(readr)
library(data.table)

#load data
dta <- read_delim("E:/Projects/Proteomics/IontrapScanRange/Reduced_human_iontrapscanrange_ITMS.txt","\t", col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
dta[dta$X1 %like% "<",]

filepath <- "E:/Projects/Proteomics/IontrapScanRange/Reduced_human_iontrapscanrange_ITMS.txt"

processFile = function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( line == "" ) {
      break
    }
    #if(length(grep("<", line)) > 0){
    #if(length(grepl("[[:digit:]]", line)) > 0){
    if(is.character(line) == TRUE){
      if(length(grep("<", line)) > 0){
        break
        } else {
          print(line)
        }
    
    } #else {TRUE}
    
  }
  
  close(con)
}


processFile = function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    #if(length(grep("<", line)) > 0){
    #if(length(grepl("[[:digit:]]", line)) > 0){
    if(is.numeric(line) == TRUE){
      print(line)
    } else {TRUE}
    
  }
  
  close(con)
}

dataList <- list()
header <- character()
add <- character()



i <- 1
for (i in 1:nrow(dta)){
  if(length(grep("<", dta[i,])) > 0){
    add[i] <- dta[i,]
    header[i] <- add[i]
  }
}




processFile(filepath)


while(length(oneLine <- readLines(fileConn, n = i, warn = FALSE)) > 0) {

  if(grepl("<", oneLine[i])){
    dataList[[i]] <- oneLine
    i <- i + 1
  }else(i <- i + 1)
}
#}






