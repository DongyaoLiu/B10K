library(ggplot2)
library(pastecs)
library(tidyverse)
library(reshape2)
library(pastecs)

group_information<-read.table("col_name.txt")
exprs<-read.table("exprs_table.txt")
group_information$species <- factor(group_information$V2)
group_information$tissues <- factor(group_information$V1)
group_information$individuals <- factor(group_information$V4)
ERCC_table <- exprs[grep("^ERCC", rownames(exprs)), ]
ERCC_table_melt <- reshape2::melt(ERCC_table)
ERCC_table_melt$ERCC_index <- factor(rep(1:92, 129))
ERCC_table_melt$ERCC_group <- factor(rep(c(rep(1,23),rep(2,23),rep(3,23),rep(4,23)),129))
for (j in c("species","tissues", "individuals")){
  tmp = list()
  for (i in group_information[ ,j]){
    tmp=append(tmp,rep(i,92),length(tmp))
  }
  ERCC_table_melt[ ,j] <- factor(unlist(tmp))
} 
###
for (i  in c(1:4)) {
  part_table<-ERCC_table_melt[grep(i,ERCC_table_melt$ERCC_group), ]
  ggplot(data = part_table, aes(x=part_table$ERCC_index, y=part_table$value)) +
    geom_point() + 
    geom_line(group=part_table$ERCC_index) +
    facet_grid(part_table$species~part_table$tissues)
  ggsave(paste(i,"pdf",sep = "."), height = 10, width = 16, units = "cm", scale = 2)
}
  
###


##############
ERCC_table.describe <- stat.desc(t(ERCC_table))
df.outlay <- c("tissues", "species","sex", "cngb_index", "ERCC_index", "count")
n=0
for (i in rownames(ERCC_table)){
n=n+1
m=0
  for (j in ERCC_table[i, ]){
    m=m+1
    b = ERCC_table.describe[11,n]
    if (j >= b){
      group_tmp <- unlist(strsplit(colnames(ERCC_table)[m], '_'))
      group_tmp <- unlist(c(group_tmp,rownames(ERCC_table)[n], ERCC_table[n,m]))
      print(group_tmp)
      df.outlay <- rbind(df.outlay,group_tmp)
    }
  }
}
rownames(df.outlay) <- NULL
df.outlay <- as.data.frame(df.outlay)
colnames(df.outlay) <- df.outlay[1, ]
df.outlay <- df.outlay[-1, ]
################

######
for (i in c(1:5)) {
  df.outlay[ , i] <- factor(df.outlay[ ,i])
}
df.outlay <- df.outlay[ ,-6]

for (i in c(1:5)) {
  tmp <- table(df.outlay[ ,i])
  plot(tmp)
}
######
avg_num_ercc = as.data.frame(table(df.outlay$cngb_index))
number <- unique(df.outlay$cngb_index)

######
df <- c("tissues", "species", "sex", "cngb")
for (i in ERCC_table_23_melt[ ,1]){
  df<-rbind(df, unlist(strsplit(i, "_")))
}
#############MASS
library(MASS)
test_df<-rep(1,2100)
scale_factor<-c(rep(1,3), rep(1.2,3), rep(0.8,3), rep(2,3), rep(5,5))
n=0
for (i in c(paste("A",c(1:3), sep = "_"), paste("B",c(1:3), sep = "_"), paste("C",c(1:3), sep = "_"), paste("D",c(1:3), sep = "_"))) {
  n=n+1
  tmp_name_1<-paste("test_gene", i, sep = "_")
  assign(tmp_name_1, 2^rnorm(2000,10,2.5)*scale_factor[n])
  tmp_name_2<-paste("test_ercc", i, sep = "_")
  assign(tmp_name_2, 2^rnorm(100,10,2.5))
  assign(i, c(get(tmp_name_1),get(tmp_name_2)))
  test_df<-cbind(test_df,get(i))
}
test_df <- test_df[ ,-1]
colnames(test_df)<-c(paste("A",c(1:3), sep = "_"), paste("B",c(1:3), sep = "_"), paste("C",c(1:3), sep = "_"), paste("D",c(1:3), sep = "_"))

#######RUV
library(RUVSeq)









 


