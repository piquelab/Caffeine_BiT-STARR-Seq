##CB 12/8/21: Script to manually do QuASAR meta-analysis
library(data.table)
library(dplyr)
library(qqman)
library(ggplot2)
water<-fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/Individual_bitstarr_water.txt") #Change as needed
a<-group_by(water, identifier) %>% summarise(m = mean(betas),sd=sd(betas),n=n(),se=sd/sqrt(n))
a$beta_z <-a$m/a$se #Make a Z score (mean/standard error)
a$new_p<-2*pnorm(-abs(a$beta_z)) #Generate p values
p<-ggplot(data=a,aes(x=new_p)) + geom_histogram()
ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/water_ase_manually_calc_pval_histo.png",p)
png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/water_ase_manually_calc_qq.png")
qq(a$new_p)
dev.off()
