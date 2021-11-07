##FastQC v0.11.9 (win/Linux file)

#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --job-name=fastq_analyse
#SBATCH --time=01:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

mkdir /data/courses/rnaseq/breastcancer_de/maelle_workspace/Quality_check/result ##create new directory that will contain the result 
cd /data/courses/rnaseq/breastcancer_de/reads/ ## enter the directory which contains the reads


for i in `ls Normal*`; do ## loop over files from the reads directory that begin with Normal => Change if you want to check another subtype of human breastcancer.
        cd /data/courses/rnaseq/breastcancer_de/maelle_workspace/Quality_check/FastQC ## Enter FastQC directory which contains the executable file for fastQC
        ./fastqc /data/courses/rnaseq/breastcancer_de/reads/$i -o ../result/ ## Launch analysis on one of the reads file and save output into the new result directory
done
