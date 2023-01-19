#!/usr/bin/env nextFlow


// Base params
runfolder_rna = params.runfolder_rna
runfolder_atac = params.runfolder_atac
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



println "============================="
println ">>> Nextflow 10x MultiOme pipeline >>>"
println ""
println "> INPUT: "
println ""
println "> run-meta-id		: $metaid "
println "> basedir		: $basedir "
println "> runfolder rna	: $runfolder_rna "
println "> runfolder atac	: $runfolder_atac "
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
println "> mgatk		: $mgatkdir "
println ""
println "============================="




// extract RNA samplesheet info
Channel
    .fromPath(sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    .tap{infoall}
    .into { crlib_ch; cragg_ch; fqc_ch; qualimap_rna_ch }

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
	val sheet into count_ready

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

mkdir -p ${aggdir}
aggcsv=${aggdir}/${projid}_libraries.csv
echo \$aggcsv
echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > \$aggcsv

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
        file "${sid}/outs/gex_possorted_bam.bam" into qualimap_rna_go
        //file "${sid}/outs/atac_possorted_bam.bam" into qualimap_atac_go
	file "${sid}/outs/atac_possorted_bam.bam" into atac_bam_ch

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
/suffolk/WorkGenomicsE/mn367/tools/cellranger-arc-2.0.1/cellranger-arc count \\
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
	"""

}




process summarize_count {

	tag "${projid}"

	input:
	val metrics from count_metrics.collect()

	output:
	val "x" into run_summarize

	"""
cd $outdir
mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger
	"""
}



process qualimap_rna_run {

	tag "${sid}-${projid}"

	input:
        val RNA_BAM from qualimap_rna_go
        set sid, projid, ref, lib, pair from qualimap_rna_ch

	output:
	val "y" into qualimap_rna_done


	"""
if [ $ref == "Human" ] || [ $ref == "human" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ]
then
        GTF=${params.custom_genome}/genes/genes.gtf.gz
else
        echo ">SPECIES NOT RECOGNIZED!"
        genome="ERR"
fi


mkdir -p  ${qcdir}/qualimap

/usr/mbu/software/qualimap/qualimap-2.2.1/qualimap rnaseq -gtf \$GTF -bam $RNA_BAM -s -outdir ${qcdir}/qualimap --java-mem-size=50G

	"""
}




/*
process qualimap_atac_run {

        tag "${sid}-${projid}"

        input:
        val ATAC_BAM from qualimap_atac_go
        set sid, projid, ref, lib, pair from qualimap_atac_ch

        output:
        val "y" into qualimap_atac_done


        """
if [ $ref == "Human" ] || [ $ref == "human" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "mouse" ] || [ $ref == "Mouse" ]
then
        GTF="/suffolk/WorkGenomicsE/mn367/Genomes/CellRanger_Genomes/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf"
elif [ $ref == "custom"  ] || [ $ref == "Custom" ]
then
        GTF=${params.custom_genome}/genes/genes.gtf.gz
else
        echo ">SPECIES NOT RECOGNIZED!"
        genome="ERR"
fi



mkdir -p  ${qcdir}/qualimap
module --ignore-cache load qualimap

qualimap bamqc -gff \$GTF -bam $ATAC_BAM -outdir ${qcdir}/qualimap --java-mem-size=50G

        """
}

*/





// aggregation
process gen_aggCSV {

	tag "${sid}_${projid}"

	input:
	set sid, projid, ref, lib, pair from cragg_ch

	output:
	set projid, ref into craggregate

	when:
	lib == 'rna'
	run_aggregate="y"

	"""


aggcsv=${aggdir}/${projid}_libraries.csv
echo \$aggcsv
echo "${sid},${aggdir}/${sid}.atac_fragments.tsv.gz,${aggdir}/${sid}.per_barcode_metrics.csv,${aggdir}/${sid}.gex_molecule_info.h5" >> \$aggcsv

	"""
}




process aggregate {

	publishDir "${outdir}/aggregate/", mode: 'move', overwrite: true
	tag "$projid"

	input:
	set projid, ref from craggregate.unique()
	val moleculeinfo from count_agg.collect()

	when:
	run_aggregate="y"

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



/suffolk/WorkGenomicsE/mn367/tools/cellranger-arc-2.0.1/cellranger-arc aggr \
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
	mgatk_run="y"

	"""
source activate mgatk

mkdir -p ${basedir}/mgatk
mkdir -p ${basedir}/mgatk/${samplename}

mgatk tenx  -bt CB -c 12 -ub UB  -i ${atac_bam} -o ${basedir}/mgatk/${samplename} -b ${barcode} --mito-genome mm10

	"""
}




process fastqc {

	tag "${sid}-${projid}"

	input:
	val run_count from count_ready
	set sid, projid, ref, lib, pair from fqc_ch

	output:
	val projid into mqc_ch

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastqc
mkdir -p ${qcdir}/fastqc/atac
mkdir -p ${qcdir}/fastqc/rna

for file in ${fqdir}/${lib}/${sid}*fastq.gz
        do
        echo \$file
        fastqc -t ${task.cpus} \$file --outdir=${qcdir}/fastqc/$lib
done

	"""
}





process multiqc_count_run {

	tag "${metaid}"

	input:
	val x from run_summarize.collect()
	val projid from mqc_ch
	val y from qualimap_rna_done

	"""

#module load multiqc
mkdir -p ${qcdir}/multiqc
~/.conda/envs/env_nf/bin/multiqc -f --outdir ${qcdir}/multiqc -n ${metaid}_sc-multiome-10x_summary_multiqc_report.html -c ${outdir}/multiqc_config.yaml ${outdir}

echo "multiqc done"

	"""
}
