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



# Convert from bed file to vcf file using PLINK:

plink --bfile all_snpsPASS10508.LDpruned --recode vcf --out all_snpsPASS10508_filterLDpruned --allow-extra-chr









