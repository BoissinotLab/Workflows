# load tidyverse package
library(tidyverse)
var_qual <- read_delim("./xenopus_genotype_subset2.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
var_qual

a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

a + theme_light() + xlim(0, 80)

var_depth <- read_delim("./xenopus_genotype_subset2.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_depth
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

summary(var_depth$mean_depth)

a + theme_light() + xlim(0, 25)

var_miss <- read_delim("./xenopus_genotype_subset2.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

summary(var_miss$fmiss)

var_freq <- read_delim("./xenopus_genotype_subset2.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)


# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

var_freq$maf

a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

a + theme_light() + xlim(0, 0.7)

summary(var_freq$maf)

ind_depth <- read_delim("./xenopus_genotype_subset2.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ind_depth
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

ind_miss  <- read_delim("./xenopus_genotype_subset2.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ind_miss
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

ind_het <- read_delim("./xenopus_genotype_subset2.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

ind_het
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

