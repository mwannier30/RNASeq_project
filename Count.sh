#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=read_count
#SBATCH --time=15:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e


module add UHTS/Analysis/subread/2.0.1;

featureCounts -a ../Mapping/reference_genome/Homo_sapiens.GRCh38.104.gtf -o ./read_count ../Mapping/outMap/*_sorted.bam

cut -f 1,7-12 ./read_count | tail -n +1 > ./cut_read_count
