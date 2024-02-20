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
cc<-inner_join(caf,wat,by="identifier")
cc$case_z<-(cc$meta_z.x-cc$meta_z.y)/sqrt(2)
cc$case_chi2 <- cc$case_z^2
ratio=median(cc$case_chi2)/0.456
cc$case_chi2_corrected <- cc$case_chi2/ratio
cc$case_p <- pchisq(cc$case_chi2_corrected,1,lower.tail=FALSE)
cc$case_padj<-p.adjust(cc$case_p)
names(cc)<-gsub(".x",".caf",names(cc))
names(cc)<-gsub(".y",".water",names(cc))

write.table(cc,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_cASE_caffeine_Update_new_meta_filt4.txt",quote=FALSE,row.names=FALSE,sep="\t")

png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/caf_case_manually_calc_qq_new_meta_filt4.png")
qq(cc$case_p)
dev.off()

