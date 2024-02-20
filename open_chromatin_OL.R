library(dplyr)
library(data.table)
load("/wsu/home/groups/piquelab/cindy/bitStarr/bitstarr_objects_dna1_dna2.RData")
ase<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/comb_ase_new_meta_filt4.txt")
case<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_cASE_caffeine_Update_new_meta_filt4.txt")
case$rsID<-gsub("_.*","",case$identifier)
ase$rsID<-gsub("_.*","",ase$identifier)
ase$padj_min<-pmin(ase$meta_padj.x,ase$meta_padj.y,na.rm=TRUE)
wa<-inner_join(ase,SNPs,by="rsID")
wa<-wa[,c(23:25,1,19,20)]
write.table(wa,file="ASE_for_chromatin_OL.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
wc<-inner_join(case,SNPs,by="rsID")
wc<-wc[,c(24:26,1,21,19,20)]
write.table(wc,file="cASE_for_chromatin_OL.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

####now run chromatin_OL.sh

asec<-fread("ASE_open_chromatin_ol.bed")
casec<-fread("cASE_open_chromatin_ol.bed")
names(asec)<-names(wa)
names(casec)<-names(wc)

sasec<-subset(asec,asec$padj_min<0.1)
scasec<-subset(casec,casec$case_padj<0.1)
length(unique(asec$identifier)) #14065
length(unique(sasec$identifier)) #330
noc<-subset(wa,!(wa$identifier %in% asec$identifier))
length(unique(noc$identifier)) #36849
snoc<-subset(noc,noc$padj_min<0.1)
length(unique(snoc$identifier)) #359
mat<-matrix(c(330,359,13735,36490),nrow=2)
#sanity check..
length(unique(wa$identifier)) #50914
sum(mat) #50914
fisher.test(mat)
#p-value < 2.2e-16
#95 percent confidence interval:
# 2.093889 2.847603
#odds ratio
#  2.442069



##Now try new p value threshold:
scasec<-subset(casec,casec$case_p<0.0215)
length(unique(casec$identifier)) #8030
length(unique(scasec$identifier)) #211
noc<-subset(wc,!(wc$identifier %in% casec$identifier))
length(unique(noc$identifier)) #15784 
snoc<-subset(noc,noc$case_p<0.0215)
length(unique(snoc$identifier)) #358
mat<-matrix(c(211,358,7819,15426),nrow=2)
#sanity check..
length(unique(wc$identifier)) #23814
sum(mat) #23814
fisher.test(mat)
#p-value = 0.08818
#95 percent confidence interval:
# 0.9740276 1.3854542
#odds ratio
 #  1.162772
