##Trimommatic v0.39

#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --job-name=fastq_analyse
#SBATCH --time=01:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip ##download Trimmomatic

cd Trimmomatic-0.39 ##enter trimmomatic folder with the executable files

## launch trimmomatich on the HER22 fastq files (paired reads)
java -jar trimmomatic-0.39.jar PE ../../../reads/HER22_R1.fastq.gz ../../../reads/HER22_R2.fastq.gz ./HER2_1trim.fastq.gz ./HER2_1unpaired.fastq.gz ./HER2_2trim.fastq.gz ./HER2_2unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
