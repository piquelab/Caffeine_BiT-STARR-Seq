library(data.table)
pwm<-fread("/nfs/rprscratch/wwwShare/carly/caffeine_huvec_bistarr_paper/Table_S2_motif_overlap_entire_lib.txt")
case<-read.table("manual_calc_DNA_Filt100_Bitstarr_cASE_caffeine_Update_new_meta_filt4.txt",header=TRUE)
ase<-read.table("comb_ase_new_meta_filt4.txt",header=TRUE)
ase$min_padj<-pmin(ase$meta_padj.x,ase$meta_padj.y)
ase$rsID<-gsub("_.*","",ase$identifier)
case$rsID<-gsub("_.*","",case$identifier)
#OL, filter, then proptest.. ASE first
detf<-inner_join(ase,pwm,by=c("rsID"="V4"))
# filter for motifs < 100 constructs in lib
te<-as.data.frame(table(pwm$V5))
tkeep<-subset(te,te$Freq>99) #Make list of motifs that have 100+ instances in lib
ndetf<- subset(detf,detf$V5 %in% tkeep$Var1)
np<-nrow(subset(ndetf,ndetf$min_padj<0.1))/nrow(ndetf)
#np=0.01075748

stf<-as.data.frame(table(subset(ndetf,ndetf$min_padj<0.1)$V5))
colnames(stf)<-c("Var1","sig_Freq")
totalde<-as.data.frame(table(detf$V5))
totalde$exp<-np
tde<-inner_join(totalde,stf,by="Var1") 
#Need to sep sig freq and total freq
fde <- sapply(tde$Var1,function(i){
a<-prop.test(tde$sig_Freq[tde$Var1==i],tde$Freq[tde$Var1==i],tde$exp[tde$Var1==i])
c(a$p.value,a$estimate,a$conf.int)
})
det<-t(fde)
rownames(det)<-tde$Var1
det<-as.data.frame(det)
colnames(det)<-c("p","estimate","ci.low","ci.high")
#Check rfs
det$tf<-rownames(det)
#subset(det,det$tf=="MA0065.2.transfac") ns
#subset(det,det$tf=="MA0155.1.transfac") enrich
#subset(det,det$tf=="MA0116.1.transfac") ns
 #subset(det,det$tf=="MA0163.1.transfac") enrich

#nrow(det) #359 tested
s<-(subset(det,p<0.05))
#nrow(subset(s,s$estimate>np)) #44 enriched
#nrow(subset(s,s$estimate<np)) #5 depleted
#First, need to match motif to name
nam<-fread("/wsu/home/groups/piquelab/Carly/JASPAR_motif_to_TF_name_final.txt",header=FALSE)
det$motif<-rownames(det)
det$motif<-gsub(".transfac","",det$motif)
det$tf<-NULL
td<-inner_join(det,nam,by=c("motif"="V1"))
write.table(td,"/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/ASE_proptest_filt100.txt")

#Now cASE (new threshold, p<0.0215), then plot
detf<-inner_join(case,pwm,by=c("rsID"="V4"))
# filter for motifs < 100 constructs in lib
te<-as.data.frame(table(pwm$V5))
tkeep<-subset(te,te$Freq>99) #Make list of motifs that have 100+ instances in lib
ndetf<- subset(detf,detf$V5 %in% tkeep$Var1)
np<-nrow(subset(ndetf,ndetf$case_p<0.0215))/nrow(ndetf)
#np=0.02799144

stf<-as.data.frame(table(subset(ndetf,ndetf$case_p<0.0215)$V5))
colnames(stf)<-c("Var1","sig_Freq")
totalde<-as.data.frame(table(detf$V5))
totalde$exp<-np
tde<-inner_join(totalde,stf,by="Var1")
#Need to sep sig freq and total freq
fde <- sapply(tde$Var1,function(i){
    a<-prop.test(tde$sig_Freq[tde$Var1==i],tde$Freq[tde$Var1==i],tde$exp[tde$Var1==i])
    c(a$p.value,a$estimate,a$conf.int)
    })
det<-t(fde)
rownames(det)<-tde$Var1
det<-as.data.frame(det)
colnames(det)<-c("p","estimate","ci.low","ci.high")
#Check rfs
det$tf<-rownames(det)
#subset(det,det$tf=="MA0065.2.transfac") #ns
#subset(det,det$tf=="MA0155.1.transfac") #enrich
#subset(det,det$tf=="MA0116.1.transfac") #ns
 #subset(det,det$tf=="MA0163.1.transfac") #ns

#nrow(det) #417 tested
s<-(subset(det,p<0.05))
#nrow(subset(s,s$estimate>np)) #18 enriched
#nrow(subset(s,s$estimate<np)) #6 depleted
#First, need to match motif to name
det$motif<-rownames(det)
det$motif<-gsub(".transfac","",det$motif)
det$tf<-NULL
td<-inner_join(det,nam,by=c("motif"="V1"))
write.table(td,"/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/cASE_proptest_filt100_new_p_threshold.txt")
