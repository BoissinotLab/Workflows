library(tidyverse)
pca <- read_table2("./xenopus.eigenvec4", col_names = FALSE)
eigenval <- scan("./xenopus.eigenval")
pca <- pca[,-1]

pca
names(pca)[1] <- "loc"
names(pca)[2] <- "ind"
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca
pve <- data.frame(PC = 1:18, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
pve

b <- ggplot(pca, aes(PC1, PC2, colour = 'loc', label = TRUE)) + geom_point(size = 3)
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b

p<-ggplot(pca,aes(x=PC1,y=PC2,color=loc))
p<-p+geom_point(size = 6)
p
p + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
