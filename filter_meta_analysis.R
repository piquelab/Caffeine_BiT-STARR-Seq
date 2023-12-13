library(data.table)
library(dplyr)
caf<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_caf_Update_new_meta.txt")
wat<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_water_Update_new_meta.txt")

caf<-subset(caf,caf$n>3)
caf$new_padj<-p.adjust(caf$meta_p,method="BH")
wat<-subset(wat,wat$n>3)
wat$new_padj<-p.adjust(wat$meta_p,method="BH")


png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/caf_ase_manually_calc_qq_new_meta_filt4.png")
qq(caf$meta_p)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/water_ase_manually_calc_qq_new_meta_filt4.png")
qq(wat$meta_p)
dev.off()

write.table(caf,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_caffeine_Update_new_meta_filt4.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(wat,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_water_Update_new_meta_filt4.txt",quote=FALSE,row.names=FALSE,sep="\t")
#Now cASE
ccaf<-inner_join(caf,wat,by="identifier")
ccaf$case_z<-ccaf$meta_z.x-ccaf$meta_z.y
ccaf$case_p<-2 * (1 - pnorm(abs(ccaf$case_z)))
ccaf$case_padj<-p.adjust(ccaf$case_p,method="BH")
names(dc)<-gsub(".x",".caf",names(dc))
names(dc)<-gsub(".y",".water",names(dc))

write.table(ccaf,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_cASE_caffeine_Update_new_meta_filt4.txt",quote=FALSE,row.names=FALSE,sep="\t")

png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/caf_case_manually_calc_qq_new_meta_filt4.png")
qq(ccaf$case_p)
dev.off()

