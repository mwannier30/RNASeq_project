##Version: HISAT2 2.1.0

#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=genome_indexing
#SBATCH --time=03:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-source/download ##download HISAT2

cd hisat2-2.1.0/ #enter the HISAT2 folder with all the executable files

./hisat2-build -p 16 ../reference_genome/*.fa ../index_files/genome_index #build the HFM index files with the reference genome as argument
