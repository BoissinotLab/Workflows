#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=23:55:00
#SBATCH --mem=118G
#SBATCH -o PCA.%J.out
#SBATCH -e PCA.%J.err

module load gencore
module load gencore_variant_detection


VCF=/scratch/da2451/Xenopus/Boissinotlab/gvcf/xenopus_genotype_filtered2.vcf.gz

# perform linkage pruning - i.e. identify prune sites

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out xenopus

# prune and create pca

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract xenopus.prune.in \
--make-bed --pca --out xenopus














