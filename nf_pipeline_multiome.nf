#!/usr/bin/env nextFlow

// Base params
basedir = params.basedir
metaid = params.metaid

// Output dirs
outdir = params.outdir
fqdir = params.fqdir
qcdir = params.qcdir
countdir = params.countdir
aggdir = params.aggdir
metadir = params.metadir
mgatkdir = params.mgatkdir

// Read and process sample sheet
sheet = file(params.sheet)

count_matrix=params.count_matrix
mgatk_run=params.mgatk_run
run_aggregate=params.run_aggregate


// software params
cellranger_arc=params.tool_cellranger_arc
qualimap=params.tool_qualimap
fastq_screen=params.tool_fastq_screen
fastqc=params.tool_fastqc
bwa=params.tool_bwa
mgatk=params.tool_mgatk
multiqc=params.tool_multiqc


println "============================="
println ">>> Nextflow 10x MultiOme pipeline >>>"
println ""
println "> INPUT: "
println ""
println "> run-meta-id		: $metaid "
println "> basedir		: $basedir "
println "> sample-sheet		: $sheet "
println ""
println ""
println " - output directory structure "
println "> outdir               : $outdir "
println "> fastq                : $fqdir "
println "> qc                   : $qcdir "
println "> count                : $countdir "
println "> aggregated           : $aggdir "
println "> metadata             : $metadir "
println "> mgatk		            : $mgatkdir "
println ""
println "============================="



// extract RNA samplesheet info
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    .tap{infoall}
    .into { crlib_ch; cragg_ch; fqc_ch; fqs_ch }

println " > Samples to process: "
println "[Sample_ID,Sample_Project,Sample_Species,Sample_Lib,pair]"



// create library file for each sample to run in parallel
process gen_libraries_csv {

    	tag "${sid}_${projid}"

	input:
	val sheet
	set sid, projid, ref, lib, pair from crlib_ch

 	output:
	set sid, projid, ref, lib, pair into count_lib_csv

	when:
	lib == 'rna'

	"""
mkdir -p $metadir
libcsv=$metadir/${projid}_${sid}_libraries.csv

# Print header
echo 'fastqs,sample,library_type' > \$libcsv

# Print RNA entry
echo '${fqdir}/rna,$sid,Gene Expression' >> \$libcsv

# Get paired ATAC sample
atacid=\$(grep ',atac,$pair' $sheet | cut -f2 -d ',')
echo "${fqdir}/atac,$sid,Chromatin Accessibility" >> \$libcsv
  	"""
}


// count RNA + ATAC
process count {

	tag "${sid}-${projid}"
	publishDir "count", mode: "copy", overwrite: true

	input:
	set sid, projid, ref, lib, pair from count_lib_csv

	output:
	path '*'
	file "${sid}/outs/filtered_feature_bc_matrix/barcodes.tsv" into barcodes_ch
	val "${qcdir}/cellranger/${sid}.summary.csv" into count_metrics
	val "${aggdir}/${sid}.molecule_info.h5" into count_agg
	val sid into mgatk_samplename_ch
	file "${sid}/outs/atac_possorted_bam.bam" into atac_bam_ch
 	path "${sid}/outs/gex_possorted_bam.bam" into rna_bam_path_ch
  	set sid, projid, ref, lib, pair into qualimap_rna_ch

	when:
	lib == 'rna'

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0"
elif [ $ref == "mouse_NUMTs_masked" ] || [ $ref == "Mouse_NUMTs_masked" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0_MT_masked/GRCm38_NUMTs_MT_5024_13715_masked_CR_arc"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ]
then
        genome=${params.custom_genome}
else
        echo ">SPECIES NOT RECOGNIZED!"
        genome="ERR"
fi


mkdir -p ${countdir}
echo $countdir

libcsv=$metadir/${projid}_${sid}_libraries.csv

$cellranger_arc count \\
	--id=$sid \\
	--libraries=\$libcsv \\
	--reference=\$genome \\
	--localmem=150 \\
	--jobmode=local \\
	--localcores=${task.cpus}


## Copy files for aggregation

# h5 file
mkdir -p $aggdir
cp ${sid}/outs/gex_molecule_info.h5 ${aggdir}/${sid}.gex_molecule_info.h5

# Atac fragments
cp ${sid}/outs/atac_fragments.tsv.gz ${aggdir}/${sid}.atac_fragments.tsv.gz
cp ${sid}/outs/atac_fragments.tsv.gz.tbi ${aggdir}/${sid}.atac_fragments.tsv.gz.tbi

# Per Barcode metrics
cp ${sid}/outs/per_barcode_metrics.csv ${aggdir}/${sid}.per_barcode_metrics.csv

## Copy metrics file for qc

# Remove if it exists
if [ -f ${qcdir}/cellranger/${sid}.summary.csv ]; then
	rm -r ${qcdir}/cellranger/${sid}.summary.csv
fi
mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger/
cp ${sid}/outs/summary.csv ${qcdir}/cellranger/${sid}.summary.csv

## Copy to delivery folder
mkdir -p ${outdir}/summaries
mkdir -p ${outdir}/summaries/web-summaries
cp ${sid}/outs/web_summary.html ${outdir}/summaries/web-summaries/${sid}.web_summary.html

gunzip -c ${sid}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ${sid}/outs/filtered_feature_bc_matrix/barcodes.tsv

mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger
	"""
}



process qualimap_rna_run {

	tag "${sid}-${projid}"

	input:
        path RNA_BAM from rna_bam_path_ch
        set sid, projid, ref, lib, pair from qualimap_rna_ch

	output:
        val projid into qualimap_rna_done

  	when:
        lib == 'rna'

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "mouse_NUMTs_masked" ] || [ $ref == "Mouse_NUMTs_masked" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0_MT_masked/GRCm38_NUMTs_MT_5024_13715_masked_CR_arc/genes/genes.gtf"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ]
then
        GTF=${params.custom_genome}/genes/genes.gtf.gz
else
        echo ">SPECIES NOT RECOGNIZED!"
        genome="ERR"
fi


mkdir -p  ${qcdir}/qualimap

$qualimap rnaseq -gtf \$GTF -bam $RNA_BAM -s -outdir ${qcdir}/qualimap/rna/${sid} --java-mem-size=50G

scp ${qcdir}/qualimap/rna/${sid}/rnaseq_qc_results.txt  ${qcdir}/qualimap/rna/${sid}_rnaseq_qc_results.txt
scp -r ${qcdir}/qualimap/rna/${sid}/raw_data_qualimapReport  ${qcdir}/qualimap/rna/${sid}_raw_data_qualimapReport
scp ${qcdir}/qualimap/rna/${sid}/qualimapReport.html  ${qcdir}/qualimap/rna/${sid}_qualimapReport.html
scp -r ${qcdir}/qualimap/rna/${sid}/css ${qcdir}/qualimap/rna/${sid}_css
scp -r ${qcdir}/qualimap/rna/${sid}/images_qualimapReport ${qcdir}/qualimap/rna/${sid}_images_qualimapReport
	"""
}



// aggregation
process gen_aggCSV {

	tag "${sid}_${projid}"

	input:
	set sid, projid, ref, lib, pair from cragg_ch

	output:
	set projid, ref into craggregate

	when:
	lib == 'rna'

	"""

mkdir -p ${aggdir}
aggcsv=${aggdir}/${projid}_libraries.csv
if [ -f \${aggcsv} ]
then
	if grep -q $sid \$aggcsv
	then
		echo ""
	else
		echo "${sid},${aggdir}/${sid}.atac_fragments.tsv.gz,${aggdir}/${sid}.per_barcode_metrics.csv,${aggdir}/${sid}.gex_molecule_info.h5" >> \$aggcsv
	fi
else
  	echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > \$aggcsv
  	echo "${sid},${aggdir}/${sid}.atac_fragments.tsv.gz,\
  	${aggdir}/${sid}.per_barcode_metrics.csv,\
  	${aggdir}/${sid}.gex_molecule_info.h5" >> \$aggcsv

fi
	"""
}




process aggregate {

	publishDir "${outdir}/aggregate/", mode: 'copy', overwrite: true
	tag "$projid"

	input:
	set projid, ref from craggregate.unique()
	val moleculeinfo from count_agg.collect()

  	output:
  	path '*'
  	file '*/web_summary.html' into cr_web_summ_mqc

	when:
	run_aggregate == 'y'

	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0"
elif [ $ref == "mouse_NUMTs_masked" ] || [ $ref == "Mouse_NUMTs_masked" ]
then
        genome="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0_MT_masked/GRCm38_NUMTs_MT_5024_13715_masked_CR_arc"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ]
then
        genome=${params.custom_genome}
else
        echo ">SPECIES NOT RECOGNIZED!"
        genome="ERR"
fi


$cellranger_arc aggr \
   --id=${projid}_agg \
   --csv=${aggdir}/${projid}_libraries.csv \
   --normalize=depth \
   --reference=\$genome

## Copy to delivery folder
cp ${projid}_agg/outs/web_summary.html ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html
#mv ${projid}_agg/outs/cloupe.cloupe ${outdir}/summaries/cloupe/${projid}_agg_cloupe.cloupe

## Remove the gex_molecule_info.h5 files that are stored in the aggregate folder (the original files are still in count-cr/../outs
rm ${aggdir}/*h5
## Remove the barcode_metrics.csv files that are stored in the aggregate folder (the original files are still in count-cr/../outs
rm ${aggdir}/*barcode_metrics.csv
## Remove the atac_fragments.csv files that are stored in the aggregate folder (the original files are still in count-cr/../outs
rm ${aggdir}/*atac_fragments.tsv.gz*
	"""
}


process mgatk {

	input:
	val samplename from mgatk_samplename_ch
	val barcode from barcodes_ch
	val atac_bam from atac_bam_ch

	when:
	mgatk_run = 'y'

	"""

mkdir -p ${basedir}/mgatk
mkdir -p ${basedir}/mgatk/${samplename}

$mgatk tenx  -bt CB -c 12 -ub UB  -i ${atac_bam} -o ${basedir}/mgatk/${samplename} -b ${barcode} --mito-genome mm10
	"""
}


process fastqc {

	tag "${sid}-${projid}"

	input:
	set sid, projid, ref, lib, pair from fqc_ch

	output:
	val projid into mqc_ch

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastqc
mkdir -p ${qcdir}/fastqc/${lib}

if [[ ${lib} = "rna" ]]
then
  $fastqc  -t ${task.cpus} ${fqdir}/${lib}/${sid}*_R*fastq.gz --outdir=${qcdir}/fastqc/$lib
elif [[ ${lib} = "atac" ]]
then
  $fastqc  -t ${task.cpus} ${fqdir}/${lib}/${sid}*_R1*fastq.gz --outdir=${qcdir}/fastqc/$lib
  $fastqc  -t ${task.cpus} ${fqdir}/${lib}/${sid}*_R3*fastq.gz --outdir=${qcdir}/fastqc/$lib
else
  echo "error, library not set"
fi
	"""
}


process fastq_screen {

	tag "${sid}-${projid}"

	input:
	set sid, projid, ref, lib, pair from fqs_ch

	output:
	val projid into mqc2_ch

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastq_screen
mkdir -p ${qcdir}/fastq_screen/atac
mkdir -p ${qcdir}/fastq_screen/rna

if [[ ${lib} = "rna" ]]
then
  $fastq_screen --aligner bwa --outdir=${qcdir}/fastq_screen/$lib ${fqdir}/${lib}/${sid}*_R*fastq.gz
elif [[ ${lib} = "atac" ]]
then
  $fastq_screen --aligner bwa --outdir=${qcdir}/fastq_screen/$lib ${fqdir}/${lib}/${sid}*_R1_*fastq.gz
  $fastq_screen --aligner bwa --outdir=${qcdir}/fastq_screen/$lib ${fqdir}/${lib}/${sid}*_R3_*fastq.gz
else
  echo "error, library not set"
fi

#rm -rf ${outdir}/qc/fastq_screen/*/*temp_subset.fastq || true
	"""
}


process multiqc_count_run {

	tag "${metaid}"

	input:
	val projid from mqc_ch
  	val projid2 from mqc2_ch
  	val projid3 from qualimap_rna_done.unique()
  	file ('*') from cr_web_summ_mqc.collect().ifEmpty([])
  	val metrics from count_metrics.collect()

	"""
mkdir -p ${qcdir}/multiqc

rm -rf ${outdir}/qc/multiqc_config.yaml || true

touch ${outdir}/qc/multiqc_config.yaml
echo "report_header_info:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Pipeline:: 'Malwina Prater'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Analysis by:: 'Malwina Prater'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Contact E-mail:: 'mn367@mrc-mbu.cam.ac.uk'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Application Type:: 'scRNA-seq + scATAC-seq'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Project Type:: '10x multiOme'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Sequencing Platform:: 'Illumina'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - Sequencing Setup:: 'standard'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "title: 'Project ID ::: $projid'" >> ${outdir}/qc/multiqc_config.yaml
echo "subtitle: 'Mitochondrial biology Unit, MRC, University of Cambridge'" >> ${outdir}/qc/multiqc_config.yaml
echo "intro_text: 'Pipeline created by Malwina Prater, MBU Bioinformatics'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "fn_ignore_files:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*lostreads*'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*I1*'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*I2*'" >> ${outdir}/qc/multiqc_config.yaml
echo "fn_ignore_dirs:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - 'work'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - 'SC_ATAC_GEX_COUNTER_CS'" >> ${outdir}/qc/multiqc_config.yaml
echo "fn_ignore_paths:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*/outs/SC_ATAC_GEX_COUNTER_CS/*'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - 'qc/atac/*R2*'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*/qc/atac/*R2*'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '*/atac/*R2*'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "fn_clean_exts:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '.gz'" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '.fastq'" >> ${outdir}/qc/multiqc_config.yaml
echo "extra_fn_clean_trim:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - '${outdir}'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "exclude_modules:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - snippy" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "custom_logo: '/suffolk/WorkGenomicsE/mn367/Logo/LOGO_MRC_MBU_Cambridge_RGB.png'" >> ${outdir}/qc/multiqc_config.yaml
echo "custom_logo_url: 'http://www.mrc-mbu.cam.ac.uk/'" >> ${outdir}/qc/multiqc_config.yaml
echo "custom_logo_title: 'MRC Mitochondrial Biology Unit, University of Cambridge'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "use_filename_as_sample_name:" >> ${outdir}/qc/multiqc_config.yaml
echo "    - qualimap" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml
echo "top_modules:" >> ${outdir}/qc/multiqc_config.yaml
echo "  - 'cellranger'" >> ${outdir}/qc/multiqc_config.yaml
echo "" >> ${outdir}/qc/multiqc_config.yaml


$multiqc -f --outdir ${qcdir}/multiqc \
	 -n ${projid}_sc-multiome-10x_summary_multiqc_report.html \
	 -c ${qcdir}/multiqc_config.yaml \
	 -d \
	 ${qcdir} \
	 ${countdir} \
	 ${fqdir} \
	 ${outdir}/summaries/web-summaries/

echo "multiqc done"
	"""
}
