# Nextflow Pipeline for 10x MultiOme processing with measurement of heteroplasmy using mgatk


## Author: Malwina Prater   <IMG SRC="Figures/nf_logo.png" width=400px><br>


--------------

## Description

Nextflow pipeline for processing 10x MultiOme datasets with CellRanger-arc and measuring mitochondrial heteroplasmy from ATAC-seq.

- Analyze 10x arc-seq (RNA + ATAC) in one pipeline.

- Starting point are fastq files



---------------

## Processing


### Input files 

1. Sample table written in `.csv`

with formatting like:

````
Sample_ID,index,Sample_Project,Sample_Species,Sample_Lib,Sample_Pair
Sample1,SINAD11,MBU_jv361_001,mouse,rna,1
Sample1,SITTH3,MBU_jv361_001,mouse,atac,1
Sample2,SINAD11,MBU_jv361_001,mouse,rna,2
Sample2,SITTH3,MBU_jv361_001,mouse,atac,2
````


2. Prepare genome.

This is already done for Mouse genome in directories (also unmasked genome for Human too):

````
/path/to/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0
/path/to/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0_MT_masked

````

For NUMTs masking of the genome `mm10.full.backlist.bed` was downloaded from `https://github.com/caleblareau/mgatk/wiki/Increasing-coverage-from-10x-processing`.

Mask genome using script: `mask_NUMTs_genome.sh`.

Once the masked fasta is ready, use script `CR_mkref.sh` to make new CellRanger-arc genome. The config file was required to proceed with making new genome.




3. Organise directory with folder `fastq` and subfolders `atac` and `rna` (`fastq.gz` files can be soft linked):

````
fastq --|---rna ---|--Sample2_S1_L001_I1_001.fastq.gz
		|		   |--Sample2_S1_L001_I2_001.fastq.gz
		|		   |--Sample2_S1_L001_R1_001.fastq.gz
		|		   |--Sample2_S1_L001_R2_001.fastq.gz
		|
		|---atac---|--Sample2_S1_L001_I1_001.fastq.gz
				   |--Sample2_S1_L001_R1_001.fastq.gz
				   |--Sample2_S1_L001_R2_001.fastq.gz
				   |--Sample2_S1_L001_R3_001.fastq.gz
````




4. Create nextflow config file using template  `nextflow.config`, and modify directories and paths to genomes and software.


5. Create config file for MultiQC `multiqc_config.yaml`, and modify project details, so it presents clearly experiment, project and collaborators.

example config file:

````
report_header_info:
    - Analysis by:: 'Malwina Prater'
    - Contact E-mail:: 'xx@xxx.ac.uk'
    - Application Type:: 'scRNA-seq + scATAC-seq'
    - Project Type:: '10x multiOme'
    - Project Owners:: 'collaborator1 & collaborator2'
    - Sequencing Platform:: 'Illumina '
    - Sequencing Setup:: 'standard'

title: "Project ID ::: MBU_ownerId_projectNo"
subtitle: "Mitochondrial biology Unit, MRC, University of Cambridge"
intro_text: "Pipeline created by Malwina Prater, MBU Bioinformatics"

#output_fn_name: MBU_ownerId_projectNo.multiqc_report.html
#data_dir_name: MBU_ownerId_projectNo.multiqc_data

fn_ignore_files:
    - '*lostreads*'

custom_logo: '/suffolk/WorkGenomicsE/mn367/Logo/LOGO_MRC_MBU_Cambridge_RGB.png'
custom_logo_url: 'http://www.mrc-mbu.cam.ac.uk/'
custom_logo_title: 'MRC Mitochondrial Biology Unit, University of Cambridge'


````


### How to run the multiome nextflow pipeline


6. Run nextflow on MBU cluster:

change directory to project dir, where fastq files are stored (or linked).

````

module load anaconda
source activate env_nf
module load qualimap

screen -S Xnf

NXF_VER=21.10.6 nextflow run nf_pipeline_multiome.nf -config nextflow.config -resume

````




## Pipeline structure:


<IMG SRC="Figures/flowchart.png" width=1000px><br>












