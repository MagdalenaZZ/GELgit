#!/bin/bash
###BSUB -J plots
#BSUB -P Analysis
#BSUB -n 4
#BSUB -q bio
#BSUB -o run_signatures_cluster.o.%J.%I
#BSUB -e run_signatures_cluster.e.%J.%I

#paramfile=signatures_sample_from_alona_filepaths_reduced
#paramfile=signatures_sample_from_ready_to_plot_11-17-17_file_locations
#paramfile=most_recent_422_cip_api_cancer_cases_filepaths
paramfile=cip_api_cancer_cases_tumour_sample_ids_plus_tumour_type_filepaths

vcf=`sed -n ${LSB_JOBINDEX}p $paramfile | awk '{print $1}'`

samplename=${vcf##*/}
samplename=${samplename%%.*}

#the below obscure looking substring grabs the second sample name from something like LP2100666-DNA_A01_LP2100665-DNA_A01, where LP2100665-DNA_A01 is the tumour sample which must match the column header in the relevant vcf otherwise there will be an error at
#  File "/home/jambrose/signatures/cancer_additional_report/CancerPlotGeneration/vcf_tools.py", line 67, in get_var_frequency
#    sample = [i for i in variant.samples if i.sample == sample_name][0]
#IndexError: list index out of range
samplename=${samplename:18:17}
echo ${samplename}

module load python/2.7.12-withsqlite

#samplename=LP2100665-DNA_A01
#samplename=LP2100666-DNA_A01_LP2100665-DNA_A01
#vcf=/genomes/by_date/2017-10-31/XP02696702/CancerLP2100665-DNA_A01_NormalLP2100666-DNA_A01/SomaticVariations/LP2100666-DNA_A01_LP2100665-DNA_A01.somatic.vcf.gz
#outdir=/home/jambrose/signatures/illumina_encore/ready_to_plot_11_17_17/${samplename}
#outdir=/home/jambrose/signatures/illumina_encore/most_recent_422_cip_api_cancer_cases/${samplename}
outdir=/home/jambrose/signatures/illumina_encore/cip_api_cancer_cases_tumour_sample_ids_plus_tumour_type/${samplename}

mkdir -p ${outdir} 

workdir=/home/jambrose/signatures/cancer_additional_report
source /home/jambrose/signatures/resources.sh

python /home/jambrose/signatures/cancer_additional_report/process_somatic_vcf_edited.py \
 --sample_name ${samplename} \
 --somatic_vcf ${vcf}  \
 --structural_variation_vcf ${workdir}/TESTDATA/TESTINPUT/TEST.somatic.SV.vcf.gz \
 --coverage_bw_file  ${workdir}/TESTDATA/TESTINPUT/TEST_TUMOR.bwtools \
 --coverage_bw_normal ${workdir}/TESTDATA/TESTINPUT/TEST_NORMAL.bwtools \
 --output_dir ${outdir} \
 --genome_ref ${REF_FASTA} \
 --genome_version grch38 \
 --conf_cellbase ${workdir}/TESTDATA/TESTINPUT/cellbase_configuration.json \
 --canonical_txs ${workdir}/TESTDATA/TESTINPUT/LIST_OF_CANONICAL_TRANSCRIPTS.tsv 
