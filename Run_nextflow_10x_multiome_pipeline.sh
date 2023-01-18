#!/bin/sh

module load anaconda
source activate env_nf
module load qualimap

screen -S Xnf

NXF_VER=21.10.6 nextflow run orig_pipeline.nf -config nextflow.config -resume
#NXF_VER=21.10.6 nextflow run orig_pipeline.nf -config nextflow.config --resume -with-dag flowchart.png

