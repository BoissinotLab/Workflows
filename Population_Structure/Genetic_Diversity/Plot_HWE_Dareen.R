################ HWE ################

HWE.all <- read.table("output_hwe_noMissing.hwe",header=T)
names(HWE.all)
attach(HWE.all)
summary(HWE.all)
levels(CHROM)

plot(HWE.all$P_HWE,xlab="Position on chr1-10LS",ylab="HWE P-value", col = ifelse(HWE.all$P_HWE < 0.05,'red','black'))

plot(HWE.all$P_HET_EXCESS,xlab="Position on chr1-10LS",ylab="HWE P-Het Excess", col = ifelse(HWE.all$P_HET_EXCESS < 0.05,'red','black'))

plot(HWE.all$P_HET_DEFICIT,xlab="Position on chr1-10LS",ylab="HWE P-Het Deficient", col = ifelse(HWE.all$P_HET_DEFICIT < 0.05,'red','black'))

plot(HWE.all$ChiSq_HWE,xlab="Position on chr1-10LS",ylab="HWE ChiSq", col = ifelse(HWE.all$ChiSq_HWE < 0.05,'red','black'))


