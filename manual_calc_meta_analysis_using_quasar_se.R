library(data.table)
library(dplyr)
library(qqman)
library(ggplot2)

water <- fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/Individual_DNA_Filt100_Bitstarr_Water_Update_CB_230313.txt")
caf <- fread("/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/Individual_DNA_Filt100_Bitstarr_Caffeine_Update_CB_230313.txt")

# Perform the meta-analysis: https://en.wikipedia.org/wiki/Inverse-variance_weighting
water_meta <- water %>%
      # Group the data by the columns 'id',
      group_by(identifier) %>%
    # Summarize the data for each group
    summarize(
          # Calculate the weighted mean of the 'beta' column
          meta_estimate = weighted.mean(betas.beta.binom, 1 / betas_se^2),
          # Calculate the standard error of the meta-analysis estimate
          meta_se = sqrt(1 / sum(1 / betas_se^2)),
          # Count the number of observations in each group
          n = n(),
          # Provide one DNA prop.
          DNA_prop = mean(DNA_prop)
        ) %>%
      # Remove the grouping applied by 'group_by'
      ungroup() %>%
      # Add columns for the meta-analysis z-score, p-value, and adjusted p-value
      mutate(
              meta_z = (meta_estimate - qlogis(DNA_prop)) / meta_se,
              meta_p = 2 * (1 - pnorm(abs(meta_z))),
              meta_padj = p.adjust(meta_p,method="BH")
            )
caf_meta <-caf %>%
      # Group the data by the columns 'id',
      group_by(identifier) %>%
    # Summarize the data for each group
    summarize(
          # Calculate the weighted mean of the 'beta' column
          meta_estimate = weighted.mean(betas.beta.binom, 1 / betas_se^2),
          # Calculate the standard error of the meta-analysis estimate
          meta_se = sqrt(1 / sum(1 / betas_se^2)),
          # Count the number of observations in each group
          n = n(),
          # Provide one DNA prop.
          DNA_prop = mean(DNA_prop)
        ) %>%
      # Remove the grouping applied by 'group_by'
      ungroup() %>%
      # Add columns for the meta-analysis z-score, p-value, and adjusted p-value
      mutate(
              meta_z = (meta_estimate - qlogis(DNA_prop)) / meta_se,
              meta_p = 2 * (1 - pnorm(abs(meta_z))),
              meta_padj = p.adjust(meta_p,method="BH")
            )


#p<-ggplot(data=water_meta,aes(x=meta_p)) + geom_histogram(bins = 100)
#ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/water_ase_manually_calc_pval_histo_new_meta.png",p)
#p<-ggplot(data=caf_meta,aes(x=meta_p)) + geom_histogram(bins = 100)
#ggsave("/nfs/rprscratch/wwwShare/carly/CindyHUVEC/caf_ase_manually_calc_pval_histo_new_meta.png",p)

png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/water_ase_manually_calc_qq_new_meta.png")
qq(water_meta$meta_p)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/CindyHUVEC/caf_ase_manually_calc_qq_new_meta.png")
qq(caf_meta$meta_p)
dev.off()


write.table(water_meta,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_water_Update_new_meta.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(caf_meta,file="/wsu/home/groups/piquelab/Carly/Cindy_HUVEC_Analysis/updated_2303/new_meta/manual_calc_DNA_Filt100_Bitstarr_caf_Update_new_meta.txt",quote=FALSE,row.names=FALSE,sep="\t")
