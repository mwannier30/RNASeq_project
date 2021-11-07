#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --job-name=fastq_analyse
#SBATCH --time=01:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

mkdir /data/courses/rnaseq/breastcancer_de/maelle_workspace/Quality_check/result
cd /data/courses/rnaseq/breastcancer_de/reads/


for i in `ls Normal*`; do
        cd /data/courses/rnaseq/breastcancer_de/maelle_workspace/Quality_check/FastQC
        ./fastqc /data/courses/rnaseq/breastcancer_de/reads/$i -o ../result/
done
