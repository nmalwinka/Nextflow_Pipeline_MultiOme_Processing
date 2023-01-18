#!/bin/sh
#SBATCH -p cluster7
#SBATCH --ntasks=6
#SBATCH --mem=64gb
#SBATCH --time=05:00:00
#SBATCH -J bwa_index
#SBATCH --export=ALL

export OMP_NUM_THREADS=8

set -e
echo pwd

module load samtools

#conda init
#conda activate
source activate base
export JAVA_HOME=~/miniconda3/bin/java
export PATH=“$JAVA_HOME/bin:$PATH”
#alias fastq_screen='/home/mn367/software/FastQ-Screen-0.14.1/fastq_screen'

echo first mask genome with FULL blacklist!

/usr/mbu/software/bedtools/bedtools-2.26.0/bin/bedtools maskfasta -fi genome.fa -bed mm10.full.backlist.bed  -fo genome_MT_5024_13715_masked.fa

echo now index genome for bwa aligner

# Make new reference genomes
bwa index genome.fa
samtools faidx genome.fa



