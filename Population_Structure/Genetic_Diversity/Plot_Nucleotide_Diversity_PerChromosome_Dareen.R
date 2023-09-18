
# Chr1S:
pi.chr1S <- subset(pi.all, CHROM == "chr1S", PI)
nrow(pi.chr1S)
summary(pi.chr1S)
plot(pi.chr1S$PI,xlab="Position on chr1S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr1S$PI < 4.191e-05,'black','red'))
plot(pi.chr1S$PI,xlab="Position on chr1S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr1S$PI < 5.618e-05,'black','blue'))

# Chr2L:
pi.chr2L <- subset(pi.all, CHROM == "chr2L", PI)
nrow(pi.chr2L)
summary(pi.chr2L)
plot(pi.chr2L$PI,xlab="Position on chr2L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr2L$PI < 4.524e-05,'black','red'))
plot(pi.chr2L$PI,xlab="Position on chr2L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr2L$PI < 7.006e-05,'black','blue'))


# Chr2S:
pi.chr2S <- subset(pi.all, CHROM == "chr2S", PI)
nrow(pi.chr2S)
summary(pi.chr2S)
plot(pi.chr2S$PI,xlab="Position on chr2S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr2S$PI < 4.444e-05,'black','red'))
plot(pi.chr2S$PI,xlab="Position on chr2S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr2S$PI < 6.764e-05,'black','blue'))

# Chr3L:
pi.chr3L <- subset(pi.all, CHROM == "chr3L", PI)
nrow(pi.chr3L)
summary(pi.chr3L)
plot(pi.chr3L$PI,xlab="Position on chr3L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr3L$PI < 4.365e-05,'black','red'))
plot(pi.chr3L$PI,xlab="Position on chr3L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr3L$PI < 6.776e-05,'black','blue'))

# Chr3S:
pi.chr3S <- subset(pi.all, CHROM == "chr3S", PI)
nrow(pi.chr3S)
summary(pi.chr3S)
plot(pi.chr3S$PI,xlab="Position on chr3S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr3S$PI < 4.191e-05,'black','red'))
plot(pi.chr3S$PI,xlab="Position on chr3S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr3S$PI < 5.465e-05,'black','blue'))

# Chr4L:
pi.chr4L <- subset(pi.all, CHROM == "chr4L", PI)
nrow(pi.chr4L)
summary(pi.chr4L)
plot(pi.chr4L$PI,xlab="Position on chr4L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr4L$PI < 4.444e-05,'black','red'))
plot(pi.chr4L$PI,xlab="Position on chr4L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr4L$PI < 6.762e-05,'black','blue'))

# Chr4S:
pi.chr4S <- subset(pi.all, CHROM == "chr4S", PI)
nrow(pi.chr4S)
summary(pi.chr4S)
plot(pi.chr4S$PI,xlab="Position on chr4S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr4S$PI < 4.333e-05,'black','red'))
plot(pi.chr4S$PI,xlab="Position on chr4S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr4S$PI < 6.547e-05,'black','blue'))

# Chr5L:
pi.chr5L <- subset(pi.all, CHROM == "chr5L", PI)
nrow(pi.chr5L)
summary(pi.chr5L)
plot(pi.chr5L$PI,xlab="Position on chr5L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr5L$PI < 4.333e-05,'black','red'))
plot(pi.chr5L$PI,xlab="Position on chr5L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr5L$PI < 6.591e-05,'black','blue'))

# Chr5S:
pi.chr5S <- subset(pi.all, CHROM == "chr5S", PI)
nrow(pi.chr5S)
summary(pi.chr5S)
plot(pi.chr5S$PI,xlab="Position on chr5S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr5S$PI < 4.333e-05,'black','red'))
plot(pi.chr5S$PI,xlab="Position on chr5S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr5S$PI < 6.591e-05,'black','blue'))

# Chr6L:
pi.chr6L <- subset(pi.all, CHROM == "chr6L", PI)
nrow(pi.chr6L)
summary(pi.chr6L)
plot(pi.chr6L$PI,xlab="Position on chr6L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr6L$PI < 4.333e-05,'black','red'))
plot(pi.chr6L$PI,xlab="Position on chr6L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr6L$PI < 6.502e-05,'black','blue'))

# Chr6S:
pi.chr6S <- subset(pi.all, CHROM == "chr6S", PI)
nrow(pi.chr6S)
summary(pi.chr6S)
plot(pi.chr6S$PI,xlab="Position on chr6S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr6S$PI < 4.191e-05,'black','red'))
plot(pi.chr6S$PI,xlab="Position on chr6S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr6S$PI < 5.764e-05,'black','blue'))

# Chr7L:
pi.chr7L <- subset(pi.all, CHROM == "chr7L", PI)
nrow(pi.chr7L)
summary(pi.chr7L)
plot(pi.chr7L$PI,xlab="Position on chr7L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr7L$PI < 4.333e-05,'black','red'))
plot(pi.chr7L$PI,xlab="Position on chr7L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr7L$PI < 6.816e-05,'black','blue'))


# Chr7S:
pi.chr7S <- subset(pi.all, CHROM == "chr7S", PI)
nrow(pi.chr7S)
summary(pi.chr7S)
plot(pi.chr7S$PI,xlab="Position on chr7S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr7S$PI < 4.333e-05,'black','red'))
plot(pi.chr7S$PI,xlab="Position on chr7S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr7S$PI < 6.372e-05,'black','blue'))

# Chr8L:
pi.chr8L <- subset(pi.all, CHROM == "chr8L", PI)
nrow(pi.chr8L)
summary(pi.chr8L)
plot(pi.chr8L$PI,xlab="Position on chr8L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr8L$PI < 4.365e-05,'black','red'))
plot(pi.chr8L$PI,xlab="Position on chr8L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr8L$PI < 6.796e-05,'black','blue'))

# Chr8S:
pi.chr8S <- subset(pi.all, CHROM == "chr8S", PI)
nrow(pi.chr8S)
summary(pi.chr8S)
plot(pi.chr8S$PI,xlab="Position on chr8S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr8S$PI < 4.333e-05,'black','red'))
plot(pi.chr8S$PI,xlab="Position on chr8S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr8S$PI < 6.287e-05,'black','blue'))

# Chr9-10L:
pi.chr9_10L <- subset(pi.all, CHROM == "chr9_10L", PI)
nrow(pi.chr9_10L)
summary(pi.chr9_10L)
plot(pi.chr9_10L$PI,xlab="Position on chr9_10L",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr9_10L$PI < 4.444e-05,'black','red'))
plot(pi.chr9_10L$PI,xlab="Position on chr9_10L",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr9_10L$PI < 6.835e-05,'black','blue'))

# Chr9-10S:
pi.chr9_10S <- subset(pi.all, CHROM == "chr9_10S", PI)
nrow(pi.chr9_10S)
summary(pi.chr9_10S)
plot(pi.chr9_10S$PI,xlab="Position on chr9_10S",ylab="Nucleotide diversity (Pi) - Above Median in Red", col = ifelse(pi.chr9_10S$PI < 4.333e-05,'black','red'))
plot(pi.chr9_10S$PI,xlab="Position on chr9_10S",ylab="Nucleotide diversity (Pi) - Above Mean in blue", col = ifelse(pi.chr9_10S$PI < 6.542e-05,'black','blue'))

