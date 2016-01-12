#merge list of tables                  
MyMerge<- function(x, y){
  df<- merge(x, y, by= 0, all.x= T, all.y= T)
  rownames(df)<- df$Row.names
  df$Row.names<- NULL
  return(df)
}

#newtable <- Reduce(MyMerge, list())

WT_3<-read.table("./WT_3_count.txt", row.names=1)
WT_4<-read.table("./WT_4_count.txt", row.names=1)
WT_5<-read.table("./WT_5_count.txt", row.names=1)
WT_6<-read.table("./WT_6_count.txt", row.names=1)
names(WT_3)<-("WT_3")
names(WT_4)<-("WT_4")
names(WT_5)<-("WT_5")
names(WT_6)<-("WT_6")
WT <-Reduce(MyMerge, list(WT_3, WT_4, WT_5, WT_6))

ref4.3_1<-read.table("./ref4-3_1_count.txt", row.names=1)
ref4.3_2<-read.table("./ref4-3_2_count.txt", row.names=1)
ref4.3_3<-read.table("./ref4-3_3_count.txt", row.names=1)
ref4.3_4<-read.table("./ref4-3_4_count.txt", row.names=1)
names(ref4.3_1)<-("ref4.3_1")
names(ref4.3_2)<-("ref4.3_2")
names(ref4.3_3)<-("ref4.3_3")
names(ref4.3_4)<-("ref4.3_4")
ref4.3 <-Reduce(MyMerge, list(ref4.3_1, ref4.3_2, ref4.3_3, ref4.3_4))

med5_1<-read.table("./med5_1_count.txt", row.names=1)
med5_2<-read.table("./med5_2_count.txt", row.names=1)
med5_3<-read.table("./med5_3_count.txt", row.names=1)
med5_4<-read.table("./med5_4_count.txt", row.names=1)
names(med5_1)<-("med5_1")
names(med5_2)<-("med5_2")
names(med5_3)<-("med5_3")
names(med5_4)<-("med5_4")
med5 <-Reduce(MyMerge, list(med5_1, med5_2, med5_3, med5_4))


med23_1<-read.table("./med23_1_count.txt", row.names=1)
med23_2<-read.table("./med23_2_count.txt", row.names=1)
med23_3<-read.table("./med23_3_count.txt", row.names=1)
med23_4<-read.table("./med23_4_count.txt", row.names=1)
names(med23_1)<-("med23_1")
names(med23_2)<-("med23_2")
names(med23_3)<-("med23_3")
names(med23_4)<-("med23_4")
med23 <-Reduce(MyMerge, list(med23_1, med23_2, med23_3, med23_4))

med2_1<-read.table("./med2_1_count.txt", row.names=1)
med2_2<-read.table("./med2_2_count.txt", row.names=1)
med2_3<-read.table("./med2_3_count.txt", row.names=1)
med2_4<-read.table("./med2_4_count.txt", row.names=1)
names(med2_1)<-("med2_1")
names(med2_2)<-("med2_2")
names(med2_3)<-("med2_3")
names(med2_4)<-("med2_4")
med2 <-Reduce(MyMerge, list(med2_1, med2_2, med2_3, med2_4))

med16_1<-read.table("./med16_1_count.txt", row.names=1)
med16_2<-read.table("./med16_2_count.txt", row.names=1)
med16_3<-read.table("./med16_3_count.txt", row.names=1)
med16_4<-read.table("./med16_4_count.txt", row.names=1)
names(med16_1)<-("med16_1")
names(med16_2)<-("med16_2")
names(med16_3)<-("med16_3")
names(med16_4)<-("med16_4")
med16 <-Reduce(MyMerge, list(med16_1, med16_2, med16_3, med16_4))

rfs33_2_1<-read.table("./rfs33_2_1_count.txt", row.names=1)
rfs33_2_2<-read.table("./rfs33_2_2_count.txt", row.names=1)
rfs33_2_3<-read.table("./rfs33_2_3_count.txt", row.names=1)
rfs33_2_4<-read.table("./rfs33_2_4_count.txt", row.names=1)
names(rfs33_2_1)<-("rfs33.2_1")
names(rfs33_2_2)<-("rfs33.2_2")
names(rfs33_2_3)<-("rfs33.2_3")
names(rfs33_2_4)<-("rfs33.2_4")
rfs33.2 <-Reduce(MyMerge, list(rfs33_2_1, rfs33_2_2, rfs33_2_3, rfs33_2_4))


##
rfs7.4_1<-read.table("./rfs7-4_med2_1_count.txt", row.names=1)
rfs7.4_2<-read.table("./rfs7-4_med2_2_count.txt", row.names=1)
rfs7.4_3<-read.table("./rfs7-4_med2_3_count.txt", row.names=1)
rfs7.4_4<-read.table("./rfs7-4_med2_4_count.txt", row.names=1)
names(rfs7.4_1)<-("rfs7.4_1")
names(rfs7.4_2)<-("rfs7.4_2")
names(rfs7.4_3)<-("rfs7.4_3")
names(rfs7.4_4)<-("rfs7.4_4")
rfs7.4 <-Reduce(MyMerge, list(rfs7.4_1, rfs7.4_2, rfs7.4_3, rfs7.4_4))


rfs2.2_1<-read.table("./rfs2-2_rfs9-4_1_count.txt", row.names=1)
rfs2.2_2<-read.table("./rfs2-2_rfs9-4_2_count.txt", row.names=1)
rfs2.2_3<-read.table("./rfs2-2_rfs9-4_3_count.txt", row.names=1)
rfs2.2_4<-read.table("./rfs2-2_rfs9-4_4_count.txt", row.names=1)
names(rfs2.2_1)<-("rfs2.2_1")
names(rfs2.2_2)<-("rfs2.2_2")
names(rfs2.2_3)<-("rfs2.2_3")
names(rfs2.2_4)<-("rfs2.2_4")
rfs2.2 <-Reduce(MyMerge, list(rfs2.2_1, rfs2.2_2, rfs2.2_3, rfs2.2_4))


rfs7.3_1<-read.table("./rfs7-3_rfs8-5_1_count.txt", row.names=1)
rfs7.3_2<-read.table("./rfs7-3_rfs8-5_2_count.txt", row.names=1)
rfs7.3_3<-read.table("./rfs7-3_rfs8-5_3_count.txt", row.names=1)
rfs7.3_4<-read.table("./rfs7-3_rfs8-5_4_count.txt", row.names=1)
names(rfs7.3_1)<-("rfs7.3_1")
names(rfs7.3_2)<-("rfs7.3_2")
names(rfs7.3_3)<-("rfs7.3_3")
names(rfs7.3_4)<-("rfs7.3_4")
rfs7.3 <-Reduce(MyMerge, list(rfs7.3_1, rfs7.3_2, rfs7.3_3, rfs7.3_4))


allcounts <-Reduce(MyMerge, list(WT, rfs33.2, med2, med5, med16, ref4.3, med23, rfs2.2, rfs7.3, rfs7.4))

justcounts <-allcounts[6:33607,]

write.csv(allcounts, "./allcounts.csv")
write.csv(justcounts, "./justcounts.csv")
