install.packages("phangorn", dependencies=TRUE)
library(phangorn)
# install latest development version needs devtools
install.packages("devtools", dependencies=TRUE)
library(usethis)
library(devtools)
library(ape)
library(phangorn)
install_github("KlausVigo/phangorn")


## RAxML best-known tree with bipartition support (from previous analysis)
raxml.tree1 <- read.tree(file.path("RAxML_bipartitions.T1"))
raxml.tree2 <- read.tree(file.path("RAxML_bestTree.T1"))
raxml.tree3 <- read.tree(file.path("RAxML_bipartitionsBranchLabels.T1"))


## RAxML bootstrap trees (from previous analysis)
raxml.bootstrap <- read.tree(file.path("RAxML_bootstrap.T1"))

par(mfrow=c(1,1), mar=c(1,1,1,1)) # Setting plot parameters for one plot
par(mfrow=c(1,3), mar=c(1,1,1,1)) # Setting plot parameters for two plots

### Plotting trees with support values:
##  RAxML
plot(raxml.tree1)
nodelabels(raxml.tree1$node.label, adj = c(1, 0), frame = "none")
plot(raxml.tree2)
nodelabels(raxml.tree2$node.label, adj = c(1, 0), frame = "none")
plot(raxml.tree3)
nodelabels(raxml.tree2$node.label, adj = c(1, 0), frame = "none")



# Show the correspondingly labelled tree and network in R
plot(raxml.tree1, "u", rotate.tree = 180, cex=.7) 
edgelabels(raxml.tree$edge[,2],col="blue", frame="none", cex=.7)
edgelabels(raxml.tree$edge[,2],col="red", frame="none", cex=.7)





