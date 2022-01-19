#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=33G
#SBATCH --job-name=genome_indexing
#SBATCH --time=20:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

module add UHTS/Aligner/hisat/2.2.1 ##import HISAT2 module
module add UHTS/Analysis/samtools/1.10 ##import samtools module

mkdir outMap ##create a new folder for the results

for reads in Normal HER2; do ##nested loop over the samples name
        for i in 1 2 3; do
                touch ./outMap/${reads}${i}_hisat2.sam ##create a new empty samfile with the sample name
                ##launch the mapping with HISAT2 with paired reads and indexes as argument and print the result in the file created above
                hisat2 -q -p 4 -x ./index_files/genome_index -1 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R1.fastq.gz -2 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R2.fastq.gz -S ./outMap/${reads}${i}_hisat2.sam
                samtools view -hbS ./outMap/${reads}${i}_hisat2.sam > ./outMap/${reads}${i}.bam ##convert sam files into bam files
                samtools sort -m 25G -@ 4 -o ./outMap/${reads}${i}_sorted.bam -T temp ./outMap/${reads}${i}.bam ##sort the resulting bam files by genomic coordinates
                samtools index ./outMap/${reads}${i}_sorted.bam ##index the resulting sorted bam files        
        done
done
