#!/bin/bash

#SBATCH -J 1_RepeatMasker
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --mem=105GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ieh211@nyu.edu

#load modules
module purge all
module load repeatmasker/4.0.7

#repeatmask
RepeatMasker -lib lacerta_bilineata-families.fa -dir . lacerta_bilineata.fa
