#!/bin/bash

#
# sets up directory structure and generates, configures and submits
# scripts to convert vcf file to hg38
# using picard LiftOverVCF
#

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
WORK_DIR=/genomes/scratch/dkasperaviciute/conversion/annotation
QUEUE=bio
THREADS=12

REFERENCE_HG19=/genomes/scratch/dkasperaviciute/conversion/input/hg19.fa
REFERENCE_HG38=/genomes/scratch/dkasperaviciute/conversion/input/hg38.fa
REFERENCE_G37=/genomes/scratch/dkasperaviciute/conversion/input/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

CHAIN_G37_HG19=/genomes/software/apps/LiftOver/GRCh37ToHg19.over.chain.gz
CHAIN_HG19_HG38=/genomes/software/apps/LiftOver/hg19ToHg38.over.chain.gz
CHAIN_HG38_HG19=/genomes/software/apps/LiftOver/hg38ToHg19.over.chain.gz
CHAIN_HG19_G37=/genomes/software/apps/LiftOver/hg19ToGRCh37.over.chain.gz

PERL_HELPER=$BASEDIR/remove_diff_backconverted.pl
NORMALISE=$BASEDIR/normalise_vcf_dedup.sh

USAGE="USAGE: qconvert2hg38.sh -i <input_dir> -n <dataset>"

while getopts "i:n:h" option; do
    case "$option" in

	i) INPUT_DIR="$OPTARG";;
	n) DATASET="$OPTARG";;
	h) echo $USAGE;;
	[?]) echo $USAGE;;

    esac
done

######
INPUT_DIR=/home/mzarowiecki/scratch/TRACERx/3_external_exomes/test_lifeover
BASEDIR=/home/mzarowiecki/git/AssemblyConversion/
DATASET=Test1
WORK_DIR=/home/mzarowiecki/scratch/TRACERx/3_external_exomes/test_work
#######



#check if required arguments are missing
if [ -z $INPUT_DIR ] || [ -z $DATASET ]
then

    echo $USAGE
    exit 1

fi



#create folders for intermediate files
RESULTS_DIR=$WORK_DIR/results/$DATASET/$TODAY
RUN_DIR=$WORK_DIR/runs/$DATASET/$TODAY
mkdir -p $RESULTS_DIR
mkdir -p $RUN_DIR

for VCF in `find $INPUT_DIR -name "*vcf.gz"`; do

	## used for UK10K as input dir contained many vcf files
#	VCF=/genomes/resources/genomeref/data/human/GRCh37_annotation/UK10K/20160215/UK10K_COHORT.20160215.sites.vcf.gz

	echo "`$NOW`"
	echo "`$NOW`###############################"
	echo "`$NOW`processing file $VCF..."
	echo "`$NOW`###############################"

	PREFIX=`basename $VCF .vcf.gz`

	SCRIPT_PATH=$RUN_DIR/convertPicard.$PREFIX.sh
	cp $BASEDIR/convertPicard.sh $SCRIPT_PATH
	chmod 770 $SCRIPT_PATH

	sed -i -e "s/#inputVCF/${VCF//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/#resultsDir/${RESULTS_DIR//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/#perlHelper/${PERL_HELPER//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/#normalise/${NORMALISE//\//\\/}/" $SCRIPT_PATH

	LOG_OUTPUT_PATH=`echo $SCRIPT_PATH | perl -pe 's/\.sh$/\.log/g'`
	echo "`${NOW}`creating and submiting job script $SCRIPT_PATH "

	bsub -q $QUEUE -R dsk -n $THREADS -o $LOG_OUTPUT_PATH $SCRIPT_PATH

done

