#!/bin/bash

#SBATCH -J repeatlandscape_divsum
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=2GB

#Load modules
module purge all
module load repeatmasker/4.0.7

#Run the calcDivergence script
../crotalus_viridis/calcDivergenceFromAlign.pl -s lacerta_bilineata.divsum lacerta_bilineata.fa.cat.gz 
