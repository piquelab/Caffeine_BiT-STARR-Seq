library(data.table)
library(dplyr)
 library(stringr)
library(ggplot2)
intact<-read.table("/nfs/rprdata/Anthony/FUNGEI/Genetics_Resub/enloc_analysis/INTACT_output/summary_file_traits_CB_20230830.txt",sep=' ',header=TRUE) #PTWAS
DF2 <-intact%>%mutate(LFDR=1-as.numeric(PCG))%>%arrange(LFDR)
x <- DF2$LFDR
FDR <- cumsum(x)/1:length(x)
DF2$FDR <- FDR
intact<-DF2


sint<-subset(intact,intact$FDR<0.1)
nrow(sint)#1121
snp<-read.table("/nfs/rprdata/Anthony/FUNGEI/Genetics_Resub/enloc_analysis/enloc_output/summary_file_SNP_CB_20231018.txt",header=TRUE)
snp$gene<-gsub(":.*","",snp$Signal)
ij<-inner_join(snp,sint,by=c("gene"="Gene"))
ij<-subset(ij,ij$SCP>0.5)
case<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/cASE_liftover_hg38.bed")
case$ID<-paste(case$V1,case$V3,sep="_")
ij$ID<-str_extract(ij$SNP, "[^_]*_[^_]*")
names(case)<-c("chr","pos","pos1","identifier","rsID","case_p","case_padj","ID")
iij<-inner_join(ij,case,by="ID")
#nrow(iij) = 1505
write.table(iij,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/INTACT_to_cASE.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
siij<-subset(iij,iij$case_padj<0.1)
# length(unique(siij$ID))
#[1] 0
siij<-subset(iij,iij$case_p<0.01)
# length(unique(siij$ID))
#[1] 8
# length(unique(siij$identifier))
#[1] 8
p <-ggplot(data=iij, aes(x=case_p)) +
        geom_histogram(position = 'identity') +
        theme_bw()
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/INTACT_to_case_pvals_quasar.png",p)
##OL with GWAS
cad<-fread("/nfs/rprdata/Anthony/FUNGEI/Genetics_Resub/zoom_plots/gwas_data/CARDIoGRAM_C4D_CAD_gwas.txt.gz")
htn<-fread("/nfs/rprdata/Anthony/FUNGEI/Genetics_Resub/zoom_plots/gwas_data/UKB_20002_1065_self_reported_hypertension_gwas.txt.gz")
names(htn)<-paste(names(htn),"HTN",sep="_")
names(cad)<-paste(names(cad),"CAD",sep="_")
 names(cad)[1]<-"ID"
 names(htn)[1]<-"ID"
 htn<-subset(htn,htn$pval_HTN<0.05)
cad<-subset(cad,cad$pval_CAD<0.05)
sgwas<-full_join(cad,htn,by="ID")
sgwas$ID<-str_extract(sgwas$ID, "[^_]*_[^_]*")
length(intersect(unique(siij$ID),unique(sgwas$ID)))
df<-inner_join(sgwas,siij,by="ID")
#Which traits?
tr<-subset(sgwas,sgwas$ID %in% df$ID)
tr[is.na(tr)] <- 1 #just to make sure it doesn't pass threshold
tr$trait<-ifelse(tr$pval_CAD<0.05,ifelse(tr$pval_HTN<0.05,"Both","CAD"),ifelse(tr$pval_HTN<0.05,"HTN","ns"))
unique(tr$trait) #Sanity check... make sure no "ns"
#[1] "Both" "HTN"
#sum(tr$trait=="Both")
#[1] 1
#sum(tr$trait=="HTN")
#[1] 7
#Now try multiple test correction within df? Updated 11/6/23 CB
df<-iij[!duplicated(iij$identifier), ]
df$new_padj<-p.adjust(df$case_p,method="BH")
subset(df,df$new_padj<0.1) #This is still the only significant variant: chr19_17057208/rs9676373_fw
