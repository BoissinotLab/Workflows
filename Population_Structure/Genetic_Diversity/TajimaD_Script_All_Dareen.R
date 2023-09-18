################ TajMD ################


taj.all <- read.table("Taj10000.Tajima.D",header=T)
names(taj.all)
attach(taj.all)
summary(taj.all)
levels(CHROM)
hist(taj.all$TajimaD,br=20)
nrow(taj.all)

