#!/bin/bash

#SBATCH -J build_database
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00

#load modules
module purge all
module load repeatmodeler/1.0.10

#build database
BuildDatabase -name "lacerta_bilineata" lacerta_bilineata.fa
