## used for merging countfiles (output from htseq-count)
## Whitney Dolan (whitneyldolan@gmail.com)
## last updated 2016-03-09

#set working directory
setwd()

#Combine with Reduce to merge lists
#Reduce(ManyMerge, list)
ManyMerge<- function(x, y){
  df<- merge(x, y, by= 0, all.x= T, all.y= T)
  rownames(df)<- df$Row.names
  df$Row.names<- NULL
  return(df)
}

#list of individual countfiles
countfiles<-list.files("./") 

#load count file with gene names as rownames, and filename as header for counts
load.file <-function(filename) {
  counttable <-read.table(filename, row.names=1, col.names=c("gene",filename))
  }

#process all count files
allcounts <-lapply(countfiles, load.file)

#merge all count files into one
mergedcounts<-Reduce(ManyMerge, allcounts)

#save merged files
write.csv(mergedcounts, "./allcounts.csv")

#save counts without summary statistics 
write.csv(mergedcounts[6:nrow(mergedcounts),], "./justcounts.csv")
