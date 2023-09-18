#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:55:00
#SBATCH --mem=118G
#SBATCH -o Genotype_stat_L.%J.out
#SBATCH -e Genotype_stat_L.%J.err

#module load gencore
#module load gencore_anaconda
#source activate /scratch/ieh211/software/tools

module load all
module load gencore
module load gencore_variant_detection

# Checking number of variants:

bcftools view -H L_genotype.vcf.gz | wc -l


# Randomly subsampling a VCF:

bcftools view L_genotype.vcf.gz | vcfrandomsample -r 0.012 > L_genotype_subset_NOfilter.vcf

# compress vcf

bcftools view L_genotype_subset_NOfilter.vcf -Oz -o L_genotype_subset_NOfilter.vcf.gz

# index vcf

bcftools index L_genotype_subset_NOfilter.vcf.gz

bcftools view -H L_genotype_subset_NOfilter.vcf.gz | wc -l

# Generating statistics from a VCF to set filtering threshold:

SET_VCF=/scratch/da2451/Xenopus/Boissinotlab/gvcf/Xenopus_laevis/All-data-pop/gvcfs/Analysis_AllData_Subgenomes/L_genotype_subset_NOfilter.vcf.gz
OUT=/scratch/da2451/Xenopus/Boissinotlab/gvcf/Xenopus_laevis/All-data-pop/gvcfs/Analysis_AllData_Subgenomes/Quality-control-stat-L 


# Calculate allele frequency:


vcftools --gzvcf $SET_VCF --freq --out $OUT  --max-alleles 2 


# Calculate mean depth per individual:

vcftools --gzvcf $SET_VCF --depth --out $OUT


# Calculate mean depth per variant:

vcftools --gzvcf $SET_VCF --site-mean-depth --out $OUT


# Calculate site quality:

vcftools --gzvcf $SET_VCF --site-quality --out $OUT


# Calculate proportion of missing data per individual:

vcftools --gzvcf $SET_VCF --missing-indv --out $OUT


# Calculate proportion of missing data per site:

vcftools --gzvcf $SET_VCF --missing-site --out $OUT


# Calculate heterozygosity:

vcftools --gzvcf $SET_VCF --het --out $OUT


# Variant per chromosome:

#vcfstats --vcf xenopus_genotype.vcf \
#        --outdir Quality_Control/ \
#        --formula 'COUNT(1) ~ CONTIG' \
#        --title 'Number of variants on each chromosome' \
#        --config Quality_Control/config.toml
