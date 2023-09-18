################ Heterozygosity ################

het.all <- read.table("output.het.het",header=T)
names(het.all)
attach(het.all)
summary(het.all)
levels(CHROM)
hist(het.all$O.HOM.,br=20)
boxplot(het.all$O.HOM.,het.all$E.HOM.,ylab="Heterozygosity")
library(ggplot2)

## Observed Homozygosity:

ggplot(het.all, aes(x = INDV, y = O.HOM.)) + geom_boxplot(aes(group= INDV)) +
  geom_jitter(position = position_jitter(0.3), aes(colour = INDV)) +
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "darkorchid4" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "mediumvioletred")

## Expected Homozygosity:

ggplot(het.all, aes(x = INDV, y = E.HOM.)) + geom_boxplot(aes(group= INDV)) +
  geom_jitter(position = position_jitter(0.3), aes(colour = INDV)) +
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "darkorchid4" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "mediumvioletred")

## F:

ggplot(het.all, aes(x = INDV, y = F)) + geom_boxplot(aes(group= INDV)) +
  geom_jitter(position = position_jitter(0.3), aes(colour = INDV)) +
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "darkorchid4" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "mediumvioletred")

## Number of sites:

ggplot(het.all, aes(x = INDV, y = N_SITES)) + geom_boxplot(aes(group= INDV)) +
  geom_jitter(position = position_jitter(0.3), aes(colour = INDV)) +
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "darkorchid4" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "mediumvioletred")
