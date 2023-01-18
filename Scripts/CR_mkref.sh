#!/bin/sh
#SBATCH -p cluster7
#SBATCH --ntasks=6
#SBATCH --mem=64gb
#SBATCH --time=05:00:00
#SBATCH -J CR_mkref
#SBATCH --export=ALL

export OMP_NUM_THREADS=8

set -e
echo pwd

source activate base
export JAVA_HOME=~/miniconda3/bin/java
export PATH=“$JAVA_HOME/bin:$PATH”

# Make new reference genomes

cellranger-arc mkref --config=/xxxx/xx/xxx/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0_MT_masked/GRCm38_NUMTs_MT_masked.config
