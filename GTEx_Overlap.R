#Need files in hg38. read in artery eQTLs overlapped with entire library:
library(data.table)
library(stringr)
e<-fread("all_GTEx_OL_entire_lib.txt")
case<-fread("cASE_liftover_hg38.bed")
ase<-fread("ASE_liftover_hg38.bed")
names(ase)<-c("chr","pos","pos1","identifier","rsID","padj_min")
names(case)<-c("chr","pos","pos1","identifier","rsID","case_p","case_padj")
e$ID<-str_extract(e$V1, "[^_]*_[^_]*")
ase$ID<-paste(ase$chr,ase$pos1,sep="_")
case$ID<-paste(case$chr,case$pos1,sep="_")
ija<-inner_join(e,ase,by="ID")
ijc<-inner_join(e,case,by="ID")
se<-subset(ija,ija$padj_min<0.1)
length(unique(ija$identifier)) #16367
length(unique(se$identifier)) #210
noc<-subset(ase,!(ase$identifier %in% ija$identifier))
length(unique(noc$identifier)) #34500
snoc<-subset(noc,noc$padj_min<0.1)
length(unique(snoc$identifier)) #470
mat<-matrix(c(210,470,16157,34030),nrow=2)
#sanity check..
length(unique(ase$identifier))
sum(mat) #23814
fisher.test(mat)
#p-value = 0.4826
#95 percent confidence interval:
# 0.7950333 1.1109028
#odds ratio
# 0.9410734

#cASE
se<-subset(ijc,ijc$case_padj<0.1)
length(unique(ijc$identifier)) #7227
length(unique(se$identifier)) #10
noc<-subset(case,!(case$identifier %in% ijc$identifier))
length(unique(noc$identifier)) #16568
snoc<-subset(noc,noc$case_padj<0.1)
length(unique(snoc$identifier)) #19
mat<-matrix(c(10,19,7217,16549),nrow=2)
#sanity check..
length(unique(case$identifier)) #23795
sum(mat) #23795
fisher.test(mat)
#p-value = 0.6866
#95 percent confidence interval:
# 0.5009984 2.7295368
#odds ratio
#  1.206834

#cASE p<0.01
se<-subset(ijc,ijc$case_p<0.01)
length(unique(ijc$identifier)) #7227
length(unique(se$identifier)) #187
noc<-subset(case,!(case$identifier %in% ijc$identifier))
length(unique(noc$identifier)) #16568
snoc<-subset(noc,noc$case_p<0.01)
length(unique(snoc$identifier)) #382
mat<-matrix(c(187,382,7040,16186),nrow=2)
#sanity check..
length(unique(case$identifier)) #23795
sum(mat) #23795
fisher.test(mat)
#p-value = 0.1964
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
 #0.9376325 1.3472573
#sample estimates:
#odds ratio
 # 1.125485

#####Also do analysis just within open chromatin regions#####
library(dplyr)
library(data.table)
library(stringr)
asec<-fread("ASE_open_chromatin_ol.bed")
casec<-fread("cASE_open_chromatin_ol.bed")
e<-fread("all_GTEx_OL_entire_lib.txt")
case<-fread("cASE_liftover_hg38.bed")
ase<-fread("ASE_liftover_hg38.bed")
names(ase)<-c("chr","pos","pos1","identifier","rsID","padj_min")
names(case)<-c("chr","pos","pos1","identifier","rsID","case_p","case_padj")
names(asec)<-c("chr","pos","pos1","identifier","rsID","padj_min")
names(casec)<-c("chr","pos","pos1","identifier","rsID","case_p","case_padj")
e$ID<-str_extract(e$V1, "[^_]*_[^_]*")
ase$ID<-paste(ase$chr,ase$pos1,sep="_")
case$ID<-paste(case$chr,case$pos1,sep="_")
ase<-subset(ase,ase$identifier %in% asec$identifier)
case<-subset(case,case$identifier %in% casec$identifier)
ija<-inner_join(e,ase,by="ID")
ijc<-inner_join(e,case,by="ID")
se<-subset(ija,ija$padj_min<0.1)
length(unique(ija$identifier)) #3461
length(unique(se$identifier)) #102
noc<-subset(ase,!(ase$identifier %in% ija$identifier))
length(unique(noc$identifier)) #10591
snoc<-subset(noc,noc$padj_min<0.1)
length(unique(snoc$identifier)) #228
mat<-matrix(c(102,228,3359,10363),nrow=2)
#sanity check..
length(unique(ase$identifier))
sum(mat) #14052
fisher.test(mat)
#p-value = 0.009564
#95 percent confidence interval:
# 1.078245 1.756716
#odds ratio
#  1.380162

#cASE
se<-subset(ijc,ijc$case_padj<0.1)
length(unique(ijc$identifier)) #7
length(unique(se$identifier)) #1976
noc<-subset(case,!(case$identifier %in% ijc$identifier))
length(unique(noc$identifier)) #6046
snoc<-subset(noc,noc$case_padj<0.1)
length(unique(snoc$identifier)) #7
mat<-matrix(c(7,7,1969,6039),nrow=2)
#sanity check..
length(unique(case$identifier)) #8022
sum(mat)
fisher.test(mat)
#p-value = 0.05435
#95 percent confidence interval:
 # 0.9166246 10.2595252
#odds ratio
 # 3.066399


#cASE p<0.01
se<-subset(ijc,ijc$case_p<0.01)
length(unique(ijc$identifier)) #1976
length(unique(se$identifier)) #59
noc<-subset(case,!(case$identifier %in% ijc$identifier))
length(unique(noc$identifier)) #6046
snoc<-subset(noc,noc$case_p<0.01)
length(unique(snoc$identifier)) #152
mat<-matrix(c(59,152,1917,5894),nrow=2)
#sanity check..
length(unique(case$identifier)) #8022
sum(mat)
fisher.test(mat)
#p-value = 0.2574
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
 #0.8641638 1.6299518
#sample estimates:
#odds ratio
 # 1.193441



