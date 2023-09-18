################ TajMD ################


# Chr1L:
Taj.chr1L <- subset(taj.all, CHROM == "chr1L", TajimaD)
nrow(Taj.chr1L)
summary(Taj.chr1L)
plot(Taj.chr1L$TajimaD,xlab="Position on chr1L",ylab="Tajima-D", col = ifelse(Taj.chr1L$TajimaD < 1.396,'black','blue'))

# Chr1S:
Taj.chr1S <- subset(taj.all, CHROM == "chr1S", TajimaD)
nrow(Taj.chr1S)
summary(Taj.chr1S)
plot(Taj.chr1S$TajimaD,xlab="Position on chr1S",ylab="Tajima-D", col = ifelse(Taj.chr1S$TajimaD < 1.396,'black','red'))

# Chr2L:
Taj.chr2L <- subset(taj.all, CHROM == "chr2L", TajimaD)
nrow(Taj.chr2L)
summary(Taj.chr2L)
plot(Taj.chr2L$TajimaD,xlab="Position on chr2L",ylab="Tajima-D", col = ifelse(Taj.chr2L$TajimaD < 1.561,'black','blue'))

# Chr2S:
Taj.chr2S <- subset(taj.all, CHROM == "chr2S", TajimaD)
nrow(Taj.chr2S)
summary(Taj.chr2S)
plot(Taj.chr2S$TajimaD,xlab="Position on chr2S",ylab="Tajima-D", col = ifelse(Taj.chr2S$TajimaD < 1.440,'black','red'))

# Chr3L:
Taj.chr3L <- subset(taj.all, CHROM == "chr3L", TajimaD)
nrow(Taj.chr3L)
summary(Taj.chr3L)
plot(Taj.chr3L$TajimaD,xlab="Position on chr3L",ylab="Tajima-D", col = ifelse(Taj.chr3L$TajimaD < 1.396,'black','blue'))

# Chr3S:
Taj.chr3S <- subset(taj.all, CHROM == "chr3S", TajimaD)
nrow(Taj.chr3S)
summary(Taj.chr3S)
plot(Taj.chr3S$TajimaD,xlab="Position on chr3S",ylab="Tajima-D", col = ifelse(Taj.chr3S$TajimaD < 1.396,'black','red'))

# Chr4L:
Taj.chr4L <- subset(taj.all, CHROM == "chr4L", TajimaD)
nrow(Taj.chr4L)
summary(Taj.chr4L)
plot(Taj.chr4L$TajimaD,xlab="Position on chr4L",ylab="Tajima-D", col = ifelse(Taj.chr4L$TajimaD < 1.494,'black','blue'))

# Chr4S:
Taj.chr4S <- subset(taj.all, CHROM == "chr4S", TajimaD)
nrow(Taj.chr4S)
summary(Taj.chr4S)
plot(Taj.chr4S$TajimaD,xlab="Position on chr4S",ylab="Tajima-D", col = ifelse(Taj.chr4S$TajimaD < 1.396,'black','red'))

# Chr5L:
Taj.chr5L <- subset(taj.all, CHROM == "chr5L", TajimaD)
nrow(Taj.chr5L)
summary(Taj.chr5L)
plot(Taj.chr5L$TajimaD,xlab="Position on chr5L",ylab="Tajima-D", col = ifelse(Taj.chr5L$TajimaD < 1.396,'black','blue'))

# Chr5S:
Taj.chr5S <- subset(taj.all, CHROM == "chr5S", TajimaD)
nrow(Taj.chr5S)
summary(Taj.chr5S)
plot(Taj.chr5S$TajimaD,xlab="Position on chr5S",ylab="Tajima-D", col = ifelse(Taj.chr5S$TajimaD < 1.494,'black','red'))

# Chr6L:
Taj.chr6L <- subset(taj.all, CHROM == "chr6L", TajimaD)
nrow(Taj.chr6L)
summary(Taj.chr6L)
plot(Taj.chr6L$TajimaD,xlab="Position on chr6L",ylab="Tajima-D", col = ifelse(Taj.chr6L$TajimaD < 1.396,'black','blue'))

# Chr6S:
Taj.chr6S <- subset(taj.all, CHROM == "chr6S", TajimaD)
nrow(Taj.chr6S)
summary(Taj.chr6S)
plot(Taj.chr6S$TajimaD,xlab="Position on chr6S",ylab="Tajima-D", col = ifelse(Taj.chr6S$TajimaD < 1.277,'black','red'))

# Chr7L:
Taj.chr7L <- subset(taj.all, CHROM == "chr7L", TajimaD)
nrow(Taj.chr7L)
summary(Taj.chr7L)
plot(Taj.chr7L$TajimaD,xlab="Position on chr7L",ylab="Tajima-D", col = ifelse(Taj.chr7L$TajimaD < 1.396,'black','blue'))

# Chr7S:
Taj.chr7S <- subset(taj.all, CHROM == "chr7S", TajimaD)
nrow(Taj.chr7S)
summary(Taj.chr7S)
plot(Taj.chr7S$TajimaD,xlab="Position on chr7S",ylab="Tajima-D", col = ifelse(Taj.chr7S$TajimaD < 1.396,'black','red'))

# Chr8L:
Taj.chr8L <- subset(taj.all, CHROM == "chr8L", TajimaD)
nrow(Taj.chr8L)
summary(Taj.chr8L)
plot(Taj.chr8L$TajimaD,xlab="Position on chr8L",ylab="Tajima-D", col = ifelse(Taj.chr8L$TajimaD < 1.396,'black','blue'))

# Chr8S:
Taj.chr8S <- subset(taj.all, CHROM == "chr8S", TajimaD)
nrow(Taj.chr8S)
summary(Taj.chr8S)
plot(Taj.chr8S$TajimaD,xlab="Position on chr8S",ylab="Tajima-D", col = ifelse(Taj.chr8S$TajimaD < 1.396,'black','red'))

# chr9_10L:
Taj.chr9_10L <- subset(taj.all, CHROM == "chr9_10L", TajimaD)
nrow(Taj.chr9_10L)
summary(Taj.chr9_10L)
plot(Taj.chr9_10L$TajimaD,xlab="Position on chr9-10L",ylab="Tajima-D", col = ifelse(Taj.chr9_10L$TajimaD < 1.498,'black','blue'))

# chr9_10S:
Taj.chr9_10S <- subset(taj.all, CHROM == "chr9_10S", TajimaD)
nrow(Taj.chr9_10S)
summary(Taj.chr9_10S)
plot(Taj.chr9_10S$TajimaD,xlab="Position on chr9-10S",ylab="Tajima-D", col = ifelse(Taj.chr9_10S$TajimaD < 1.396,'black','red'))

