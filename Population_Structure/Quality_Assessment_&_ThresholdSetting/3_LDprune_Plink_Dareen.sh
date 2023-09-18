#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=23:55:00
#SBATCH --mem=118G
#SBATCH -o LDprune_PCA_errbii.%J.out
#SBATCH -e LDprune_PCA_errbii.%J.err


module load all
module load gencore
module load gencore_variant_detection


#VCF=//scratch/da2451/Xenopus/Boissinotlab/gvcf/Xenopus_laevis/All-data-pop/gvcfs/Analysis_AllData/SambaR/LDPrune/xenopus_genotypeALL_filtered.vcf.gz

## Convert vcf to plink format wusing vcfTools and PLINK:


# generate chromosomes map:
grep -v "^#" xenopus_genotypeALL_filtered.vcf|cut -f1 | uniq| sort -V|awk '{print $0"\t"$0}' > chrom.map

# vcf >> plink:

vcftools --vcf xenopus_genotypeALL_filtered.vcf --plink --chrom-map chrom.map --out all_snpsPASS

# generate bed file:

plink --file all_snpsPASS --make-bed --aec --out all_snpsPASS


#### Linkage Disequilibrium (LD) pruning using PLINK:

# first generate a list of position to keep/remove - R^2 of 0.8 is considered "tagging":

plink --bfile all_snpsPASS --aec --indep-pairwise 10 5 0.8 --out all_snpsPASS10508.prune


# outputs two files: all_snpsPASS.prune.in (to keep) and all_snpsPASS.prune.out (to remove)

# then

plink --bfile all_snpsPASS --exclude all_snpsPASS10508.prune.prune.out --aec --out all_snpsPASS10508.LDpruned --make-bed


# Convert from bed file to vcf file using PLINK:

#plink --bfile my_file_prefix --recode vcf --out my_output_file_prefix












