#!/bin/bash
#SBATCH -J Pop_stat2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=105GB
#SBATCH -o Pop.%J.out
#SBATCH -e Pop.%J.err

module purge all
module load vcftools/intel/0.1.14
module load bcftools/intel/1.9

# Zip File:

#bcftools view xenopus_genotype_filtered3.vcf -Oz -o xenopus_genotype_filtered3.vcf.gz

#bcftools view xenopus_genotype_filtered3.vcf.idx -Oz -o xenopus_genotype_filtered3.vcf.idx.gz


# Output a Hardy-Weinberg p-value for every site in the bcf file  that  does  not  have  any missing genotypes:

vcftools --vcf xenopus_genotype_filtered3.vcf --hardy --max-missing 1.0 --out output_hwe_noMissing
vcftools --vcf xenopus_genotype_filtered3.vcf –hwe 0.001 --recode

# Output nucleotide diversity:

vcftools --vcf xenopus_genotype_filtered3.vcf --window-pi  10000 --out xenopus_NDiv

# LD:

vcftools --vcf xenopus_genotype_filtered3.vcf --hap-r2 --out output_LD

# Fst estimate from Weir and Cockerham's 1984:

vcftools --vcf xenopus_genotype_filtered3.vcf --weir-fst-pop North.txt --weir-fst-pop Middle.txt --weir-fst-pop South.txt --out popNMS_FST

# Heterozygosity: Calculates a measure of heterozygosity on a  per-individual  basis.  Specfically,  the inbreeding coefficient, F, is estimated for each individual using a method of moments.The resulting file has the suffix ".het".

vcftools --het xenopus_genotype_filtered3.vcf --out output_hetero

vcftools --vcf xenopus_genotype_filtered3.vcf --het --out output.het

# Tajm_D:

#vk tajima 5000 3000 xenopus_genotype_filtered3.vcf > output.TajimaD

vcftools --vcf xenopus_genotype_filtered3.vcf --TajimaD 10000 --out Taj10000


