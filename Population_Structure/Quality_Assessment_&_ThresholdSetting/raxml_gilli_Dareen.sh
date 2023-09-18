#!/bin/bash

#SBATCH -J raxml
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB
#SBATCH --time=240:00:00
#SBATCH -o Raxml_gilli.%J.out
#SBATCH -e Raxml_gilli.%J.err

#load modules
module load gencore
module load gencore_anaconda/3-4.0.0
source activate /home/da2451/.conda/envs/raxML

raxmlHPC -f a -m GTRGAMMA -p 12355 -x 12355 -# 100 -s xenopus_genotype_subset_gilli.min4.phy -n T1


