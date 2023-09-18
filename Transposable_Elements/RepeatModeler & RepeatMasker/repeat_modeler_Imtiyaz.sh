#!/bin/bash

#SBATCH -J repeat_modeler
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=144:00:00
#SBATCH --mem=105GB

#load modules
module purge all
module load repeatmodeler/1.0.10

#repeatmodeler
RepeatModeler -pa 24 -database lacerta_bilineata 
