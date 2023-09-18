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
module load gencore_anaconda

VCF=/scratch/da2451/Xenopus/Boissinotlab/gvcf/Xenopus_laevis/gvcf2/gvcfs/Quality_Control/Filters/xenopus_genotype_filtered.vcf.gz


# perform linkage pruning - i.e. identify prune sites

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out xenopus_pca

# prune and create pca

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract xenopus_pca.prune.in \
--make-bed --pca --out xenopus_pca














