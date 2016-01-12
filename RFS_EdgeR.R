#source("http://bioconductor.org/biocLite.R")
#biocLite()

#install edgeR package
#biocLite("edgeR")
#library(edgeR) #loads package into current session

#merge list of tables                  
MyMerge<- function(x, y){
  df<- merge(x, y, by= 0, all.x= T, all.y= T)
  rownames(df)<- df$Row.names
  df$Row.names<- NULL
  return(df)
}

#read counts file
#updated justcounts files on 10/24/15 to reflect that samples med2_4 and rfs7-4_4 were switched. 
x <- read.csv("./count files/justcounts.csv", row.names=1, stringsAsFactors=FALSE)

#put counts and other information into a DGElist object
#updated 10/24/15 to reflect that samples med2_4 and rfs7-4_4 were switched. 
group <- factor(c("WT","WT","WT","WT","rfs33.2","rfs33.2","rfs33.2","rfs33.2","med2","med2","med2","med2","med5","med5","med5","med5","med16","med16","med16","med16","ref4.3","ref4.3","ref4.3","ref4.3","med23","med23","med23","med23","rfs2.2","rfs2.2","rfs2.2","rfs2.2","rfs7.3","rfs7.3","rfs7.3","rfs7.3","rfs7.4","rfs7.4","rfs7.4","rfs7.4"))

dge <- DGEList(counts=x, group=group)
dim(dge) #33602


### Filtering the data ###

#filter out lowly expressed tags 
keep <- rowSums(cpm(dge)>1) >= 4  #keep rows where cpm is greater than 1 in at least four samples
dge <- dge[keep,]
dim(dge)   #18809
#re-compute the library size
dge$samples$lib.size <- colSums(dge$counts)


### Normalizing ###

#Compute effective library sizes using TMM normalization:
dge <- calcNormFactors(dge)
dge$samples  


### Data Exploration ### 
#An MDS plots shows distances, in terms of biological coefficient of variation (BCV), between samples:
plotMDS(dge) 

##glm model##
#construct design matrix
design <- model.matrix(~group)
#estimate dispersion
dge<- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge,design)
dge <- estimateGLMTagwiseDisp(dge,design)

#ANOVA 
fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=2:10)
topTags(lrt)

ANOVAresults<-lrt$table
adj<-p.adjust(ANOVAresults$PValue, method="BH")
ANOVAresults<-cbind(ANOVAresults,adj)


write.csv(ANOVAresults, "./ANOVAresults.csv")

##pairwise comparisons with WT##
dge <- estimateCommonDisp(dge, verbose=TRUE)
dge <- estimateTagwiseDisp(dge)

WT_33.2.et <-exactTest(dge, pair=c("WT", "rfs33.2"))
WT_ref4.3.et <-exactTest(dge, pair=c("WT", "ref4.3"))
WT_med2.et <-exactTest(dge, pair=c("WT", "med2"))
WT_med5.et <-exactTest(dge, pair=c("WT", "med5"))
WT_med16.et <-exactTest(dge, pair=c("WT", "med16"))
WT_med23.et <-exactTest(dge, pair=c("WT", "med23"))
WT_7.3.et <-exactTest(dge, pair=c("WT", "rfs7.3"))
WT_7.4.et <-exactTest(dge, pair=c("WT", "rfs7.4"))
WT_2.2.et <-exactTest(dge, pair=c("WT", "rfs2.2"))

WT_ref4.3 <-WT_ref4.3.et$table
adjP<-p.adjust(WT_ref4.3$PValue, method="BH")
WT_ref4.3<-cbind(WT_ref4.3,adjP)
colnames(WT_ref4.3)[1:4]<- paste("ref4.3", colnames(WT_ref4.3)[1:4], sep = "_")

WT_med2 <-WT_med2.et$table
adjP<-p.adjust(WT_med2$PValue, method="BH")
WT_med2<-cbind(WT_med2,adjP)
colnames(WT_med2)[1:4]<- paste("med2", colnames(WT_med2)[1:4], sep = "_")

WT_med5 <-WT_med5.et$table
adjP<-p.adjust(WT_med5$PValue, method="BH")
WT_med5<-cbind(WT_med5,adjP)
colnames(WT_med5)[1:4]<- paste("med5", colnames(WT_med5)[1:4], sep = "_")

WT_med16 <-WT_med16.et$table
adjP<-p.adjust(WT_med16$PValue, method="BH")
WT_med16<-cbind(WT_med16,adjP)
colnames(WT_med16)[1:4]<- paste("med16", colnames(WT_med16)[1:4], sep = "_")

WT_med23 <-WT_med23.et$table
adjP<-p.adjust(WT_med23$PValue, method="BH")
WT_med23<-cbind(WT_med23,adjP)
colnames(WT_med23)[1:4]<- paste("med23", colnames(WT_med23)[1:4], sep = "_")

WT_33.2 <-WT_33.2.et$table
adjP<-p.adjust(WT_33.2$PValue, method="BH")
WT_33.2<-cbind(WT_33.2,adjP)
colnames(WT_33.2)[1:4]<- paste("33.2", colnames(WT_33.2)[1:4], sep = "_")

WT_7.3 <-WT_7.3.et$table
adjP<-p.adjust(WT_7.3$PValue, method="BH")
WT_7.3<-cbind(WT_7.3,adjP)
colnames(WT_7.3)[1:4]<- paste("7.3", colnames(WT_7.3)[1:4], sep = "_")

WT_7.4 <-WT_7.4.et$table
adjP<-p.adjust(WT_7.4$PValue, method="BH")
WT_7.4<-cbind(WT_7.4,adjP)
colnames(WT_7.4)[1:4]<- paste("7.4", colnames(WT_7.4)[1:4], sep = "_")

WT_2.2 <-WT_2.2.et$table
adjP<-p.adjust(WT_2.2$PValue, method="BH")
WT_2.2<-cbind(WT_2.2,adjP)
colnames(WT_2.2)[1:4]<- paste("2.2", colnames(WT_2.2)[1:4], sep = "_")



#newtable <- Reduce(MyMerge, list())

allpairwise <-Reduce(MyMerge, list(WT_ref4.3, WT_33.2, WT_med2, WT_med5, WT_med16, WT_med23, WT_7.3, WT_7.4, WT_2.2))

write.csv(allpairwise, "./pairwise_WT.csv")

##pairwise comparisons between ref4-3 and rfs##

dge <- estimateCommonDisp(dge, verbose=TRUE)
dge <- estimateTagwiseDisp(dge)

ref4.3_33.2.et <-exactTest(dge, pair=c("ref4.3", "rfs33.2"))
ref4.3_7.3.et <-exactTest(dge, pair=c("ref4.3", "rfs7.3"))
ref4.3_7.4.et <-exactTest(dge, pair=c("ref4.3", "rfs7.4"))
ref4.3_2.2.et <-exactTest(dge, pair=c("ref4.3", "rfs2.2"))
ref4.3_med2.et <-exactTest(dge, pair=c("ref4.3", "med2"))
ref4.3_med5.et <-exactTest(dge, pair=c("ref4.3", "med5"))
ref4.3_med16.et <-exactTest(dge, pair=c("ref4.3", "med16"))
ref4.3_med23.et <-exactTest(dge, pair=c("ref4.3", "med23"))


ref4.3_33.2 <-ref4.3_33.2.et$table
adjP<-p.adjust(ref4.3_33.2$PValue, method="BH")
ref4.3_33.2<-cbind(ref4.3_33.2,adjP)
colnames(ref4.3_33.2)[1:4]<- paste("rfs33.2", colnames(ref4.3_33.2)[1:4], sep = "_")

ref4.3_rfs7.3 <-ref4.3_7.3.et$table
adjP<-p.adjust(ref4.3_rfs7.3$PValue, method="BH")
ref4.3_rfs7.3<-cbind(ref4.3_rfs7.3,adjP)
colnames(ref4.3_rfs7.3)[1:4]<- paste("rfs7.3", colnames(ref4.3_rfs7.3)[1:4], sep = "_")

ref4.3_rfs7.4 <-ref4.3_7.4.et$table
adjP<-p.adjust(ref4.3_rfs7.4 $PValue, method="BH")
ref4.3_rfs7.4 <-cbind(ref4.3_rfs7.4, adjP)
colnames(ref4.3_rfs7.4)[1:4]<- paste("rfs7.4", colnames(ref4.3_rfs7.4 )[1:4], sep = "_")

ref4.3_rfs2.2 <-ref4.3_2.2.et$table
adjP<-p.adjust(ref4.3_rfs2.2 $PValue, method="BH")
ref4.3_rfs2.2 <-cbind(ref4.3_rfs2.2, adjP)
colnames(ref4.3_rfs2.2)[1:4]<- paste("rfs2.2", colnames(ref4.3_rfs2.2)[1:4], sep = "_")

ref4.3_med2 <-ref4.3_2.2.et$table
adjP<-p.adjust(ref4.3_med2 $PValue, method="BH")
ref4.3_med2 <-cbind(ref4.3_med2, adjP)
colnames(ref4.3_med2)[1:4]<- paste("med2", colnames(ref4.3_med2)[1:4], sep = "_")

ref4.3_med5 <-ref4.3_2.2.et$table
adjP<-p.adjust(ref4.3_med5 $PValue, method="BH")
ref4.3_med5 <-cbind(ref4.3_med5, adjP)
colnames(ref4.3_med5)[1:4]<- paste("med5", colnames(ref4.3_med5)[1:4], sep = "_")

ref4.3_med16 <-ref4.3_2.2.et$table
adjP<-p.adjust(ref4.3_med16 $PValue, method="BH")
ref4.3_med16 <-cbind(ref4.3_med16, adjP)
colnames(ref4.3_med16)[1:4]<- paste("med16", colnames(ref4.3_med16)[1:4], sep = "_")

ref4.3_med23 <-ref4.3_2.2.et$table
adjP<-p.adjust(ref4.3_med23 $PValue, method="BH")
ref4.3_med23 <-cbind(ref4.3_med23, adjP)
colnames(ref4.3_med23)[1:4]<- paste("med23", colnames(ref4.3_med23)[1:4], sep = "_")


ref4.3_pairwise <-Reduce(MyMerge, list(ref4.3_33.2, ref4.3_rfs7.3, ref4.3_rfs7.4, ref4.3_rfs2.2, ref4.3_med2, ref4.3_med5, ref4.3_med16, ref4.3_med23))

write.csv(ref4.3_pairwise, "./ref4-3_pairwise.csv")

##########
TAIR10fun<-(read.csv("./references/TAIR10_functional_descriptions2.csv"))
TAIR10fun<-()
TAIR10fun<-data.frame(TAIR10fun, row.names=1)

allpairwise2 <-Reduce(MyMerge, list(allpairwise, TAIR10fun))

####UNDER CONSTRUCTION####
#Plot the log-fold-changes, highlighting the DE genes:
detags <- rownames(ydge)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#export all DEG at 5% FDR

DEG <- topTags(et, n=countDGE) 

write.csv(DEG, "/Users/Whitney/Documents/FLC/DEG_0.05FDR.csv")


