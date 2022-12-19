library(DESeq2)
library(EnhancedVolcano)
library(qqman)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(data.table)

sample<-fread("covariate_table.txt")
#Read in files
SS1<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HC1.q20.query.txt")
SS2<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HC2.q20.query.txt")
SS3<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HC3.q20.query.txt")
SS4<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HW2.q20.query.txt")
SS5<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HC2.q20.query.txt")
SS6<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HC3.q20.query.txt")

WW1<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HW1.q20.query.txt")
WW2<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HW2.q20.query.txt")
WW3<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HW3.q20.query.txt")
WW4<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HW1.q20.query.txt")
 WW5<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST5HC1.q20.query.txt")
WW6<-fread("/wsu/home/groups/piquelab/cindy/STARR/HUVECs/working/dedup.BST11HW3.q20.query.txt")

SS1$V10<-paste(SS1$V1,SS1$V2,sep=":")
SS2$V10<-paste(SS2$V1,SS2$V2,sep=":")
SS3$V10<-paste(SS3$V1,SS3$V2,sep=":")
SS4$V10<-paste(SS4$V1,SS4$V2,sep=":")
SS5$V10<-paste(SS5$V1,SS5$V2,sep=":")
SS6$V10<-paste(SS6$V1,SS6$V2,sep=":")
WW1$V10<-paste(WW1$V1,WW1$V2,sep=":")
WW2$V10<-paste(WW2$V1,WW2$V2,sep=":")
WW3$V10<-paste(WW3$V1,WW3$V2,sep=":")
WW4$V10<-paste(WW4$V1,WW4$V2,sep=":")
WW5$V10<-paste(WW5$V1,WW5$V2,sep=":")
WW6$V10<-paste(WW6$V1,WW6$V2,sep=":")


#Join all to create 1 df
SS12<-full_join(SS1,SS2)
SS13<-full_join(SS12,SS3)
SS14<-full_join(SS13,SS4)
SS15<-full_join(SS14,SS5)
SS16<-full_join(SS15,SS6)

WW12<-full_join(WW1,WW2)
WW13<-full_join(WW12,WW3)
WW14<-full_join(WW13,WW4)
WW15<-full_join(WW14,WW5)
WW16<-full_join(WW15,WW6)
#filter to remove rows where all NA from both dfs
WW16 %>% filter_all(any_vars(!is.na(.)))
SS16 %>% filter_all(any_vars(!is.na(.)))
 ID<-intersect(SS16$V10,WW16$V10)
 ID<-unique(ID)
 test<-as.data.frame(ID)
#sum fwd
SS1$V11<-SS1$V6 + SS1$V8
SS2$V11<-SS2$V6 + SS2$V8
SS3$V11<-SS3$V6 + SS3$V8
SS4$V11<-SS4$V6 + SS4$V8
SS5$V11<-SS5$V6 + SS5$V8
SS6$V11<-SS6$V6 + SS6$V8
WW1$V11<- WW1$V6 + WW1$V8
WW2$V11<- WW2$V6 + WW2$V8
WW3$V11<- WW3$V6 + WW3$V8
WW4$V11<- WW4$V6 + WW4$V8
WW5$V11<- WW5$V6 + WW5$V8
WW6$V11<- WW6$V6 + WW6$V8
#sum rev
SS1$V12<-SS1$V7 + SS1$V9
SS2$V12<-SS2$V7 + SS2$V9
SS3$V12<-SS3$V7 + SS3$V9
SS4$V12<-SS4$V7 + SS4$V9
SS5$V12<-SS5$V7 + SS5$V9
SS6$V12<-SS6$V7 + SS6$V9
WW1$V12<- WW1$V7 + WW1$V9
WW2$V12<- WW2$V7 + WW2$V9
WW3$V12<- WW3$V7 + WW3$V9
WW4$V12<- WW4$V7 + WW4$V9
WW5$V12<- WW5$V7 + WW5$V9
WW6$V12<- WW6$V7 + WW6$V9
df1<-left_join(test, SS1[,10:11],by=c("ID"="V10"))
df2<-left_join(test, SS2[,10:11],by=c("ID"="V10"))
df3<-left_join(test, SS3[,10:11],by=c("ID"="V10"))
df4<-left_join(test, SS4[,10:11],by=c("ID"="V10"))
df5<-left_join(test, SS5[,10:11],by=c("ID"="V10"))
df6<-left_join(test, SS6[,10:11],by=c("ID"="V10"))
df7<-left_join(test, WW1[,10:11],by=c("ID"="V10"))
df8<-left_join(test, WW2[,10:11],by=c("ID"="V10"))
df9<-left_join(test, WW3[,10:11],by=c("ID"="V10"))
df10<-left_join(test, WW4[,10:11],by=c("ID"="V10"))
df11<-left_join(test, WW5[,10:11],by=c("ID"="V10"))
df12<-left_join(test, WW6[,10:11],by=c("ID"="V10"))
colnames(df1)<-c("ID","S1")
colnames(df2)<-c("ID","S2")
colnames(df3)<-c("ID","S3")
colnames(df4)<-c("ID","S4")
colnames(df5)<-c("ID","S5")
colnames(df6)<-c("ID","S6")
colnames(df7)<-c("ID","W1")
colnames(df8)<-c("ID","W2")
colnames(df9)<-c("ID","W3")
 colnames(df10)<-c("ID","W4")
colnames(df11)<-c("ID","W5")
colnames(df12)<-c("ID","W6")
cf1<-left_join(test, SS1 %>% dplyr::select(V10,V12) ,by=c("ID"="V10"))
cf2<-left_join(test, SS2 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf3<-left_join(test, SS3 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf4<-left_join(test, SS4 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf5<-left_join(test, SS5 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf6<-left_join(test, SS6 %>% dplyr::select(V10,V12),by=c("ID"="V10")) 
cf7<-left_join(test, WW1 %>% dplyr::select(V10,V12) ,by=c("ID"="V10"))
cf8<-left_join(test, WW2 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf9<-left_join(test, WW3 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf10<-left_join(test, WW4 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf11<-left_join(test, WW5 %>% dplyr::select(V10,V12),by=c("ID"="V10"))
cf12<-left_join(test, WW6 %>% dplyr::select(V10,V12),by=c("ID"="V10")) 
colnames(cf1)<-c("ID","S1")
colnames(cf2)<-c("ID","S2")
colnames(cf3)<-c("ID","S3")
colnames(cf4)<-c("ID","S4")
colnames(cf5)<-c("ID","S5")
colnames(cf6)<-c("ID","S6")
colnames(cf7)<-c("ID","W1")
colnames(cf8)<-c("ID","W2")
colnames(cf9)<-c("ID","W3")
 colnames(cf10)<-c("ID","W4")
colnames(cf11)<-c("ID","W5")
colnames(cf12)<-c("ID","W6")

#Now inner_join() all
table_r<-inner_join(df1,df2)
table_r<-inner_join(table_r,df3)
table_r<-inner_join(table_r,df4)
table_r<-inner_join(table_r,df5)
table_r<-inner_join(table_r,df6)
table_r<-inner_join(table_r,df7)
table_r<-inner_join(table_r,df8)
table_r<-inner_join(table_r,df9)
table_r<-inner_join(table_r,df10)
table_r<-inner_join(table_r,df11)
table_r<-inner_join(table_r,df12)
table_r = table_r[!duplicated(table_r$ID),]
rownames(table_r)<-table_r[,1]
table_r<-table_r[-c(1)]
table_r[is.na(table_r)] <- 0
rowsToKeep = rowSums(table_r>0)>=(0.5*ncol(table_r))
table_r <- table_r[rowsToKeep, ]
 
table_a<-inner_join(cf1,cf2)
table_a<-inner_join(table_a,cf3)
table_a<-inner_join(table_a,cf4)
table_a<-inner_join(table_a,cf5)
table_a<-inner_join(table_a,cf6)
table_a<-inner_join(table_a,cf7)
table_a<-inner_join(table_a,cf8)
table_a<-inner_join(table_a,cf9)
table_a<-inner_join(table_a,cf10)
table_a<-inner_join(table_a,cf11)
table_a<-inner_join(table_a,cf12)
table_a = table_a[!duplicated(table_a$ID),]
rownames(table_a)<-table_a[,1]
table_a<-table_a[-c(1)]
table_a[is.na(table_a)] <- 0
rowsToKeep = rowSums(table_a>0)>=(0.5*ncol(table_a))
table_a <- table_a[rowsToKeep, ]
data<-DESeqDataSetFromMatrix(table_r, sample, ~allele+condition)
dataDE<-DESeq(data,)
 dres<-results(dataDE)
 summary(dres)
 
##Make df
res<-as.data.frame(dres)
png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_rev_Caffeine_qq.png")
q<-qq(res$pvalue)
dev.off()
a<-ggplot(data=res,aes(x =pvalue)) + geom_histogram()   
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_rev_Caffeine_pval_histo.png",a)
 
vsd<-varianceStabilizingTransformation(dataDE,blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
b<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
 ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_rev_Caffeine_heatmap.png",b)
 
c<-plotPCA(vsd, intgroup=c("condition", "sample"),returnData=TRUE)
percentVar <- round(100 * attr(c, "percentVar"))
d<-ggplot(c,aes(PC1,PC2,color=sample,shape=condition))+geom_point(size=5)+  xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_rev_Caffeine_pca.png",d)
e<-EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    ylab = bquote(~-Log[10] ~ italic(Adj_P)),
   pCutoff = 0.1,    FCcutoff = 0,
)
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_rev_Caffeine_volcano.png",e)
resOrdered <- res[order(res$padj),]
write.table(resOrdered,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/DESeq_Caffeine_rev.txt",quote=FALSE)

#Rerun for table_a 
data<-DESeqDataSetFromMatrix(table_a, sample, ~allele+condition)
dataDE<-DESeq(data,)
 dres<-results(dataDE)
 summary(dres)
 
##Make df
res<-as.data.frame(dres)
png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_fwd_Caffeine_qq.png")
q<-qq(res$pvalue)
dev.off()
a<-ggplot(data=res,aes(x =pvalue)) + geom_histogram()   
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_fwd_Caffeine_pval_histo.png",a)
 
vsd<-varianceStabilizingTransformation(dataDE,blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(fwd(brewer.pal(9, "Blues")) )(255)
b<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
 ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_fwd_Caffeine_heatmap.png",b)
 
c<-plotPCA(vsd, intgroup=c("condition", "sample"),returnData=TRUE)
percentVar <- round(100 * attr(c, "percentVar"))
d<-ggplot(c,aes(PC1,PC2,color=sample,shape=condition))+geom_point(size=5)+  xlab(paste0("PC1: ",percentVar[1],"% variance")) +ylab(paste0("PC2: ",percentVar[2],"% variance"))
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_fwd_Caffeine_pca.png",d)
e<-EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    ylab = bquote(~-Log[10] ~ italic(Adj_P)),
   pCutoff = 0.1,    FCcutoff = 0,
)
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/DESeq_Caffeine/HUVEC_fwd_Caffeine_volcano.png",e)
resOrdered <- res[order(res$padj),]
write.table(resOrdered,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/DESeq_Caffeine_fwd.txt",quote=FALSE)
