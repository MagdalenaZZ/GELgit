#!/bin/bash

#template of the job for vcf file convertion from GRCh37 to GRCh38
#shopt -s expand_aliases
#load modules
#module load java/jdk1.8.0_45
#export JAVA_OPTS='-XX:ParallelGCThreads=3'
module load picard/2.12.1
module load bcftools/1.5
module load vt/ee9a751
NOW="date +%Y-%m-%d%t%T%t"

#number of sorting threads for normalisation
N_THREADS=12

#reference and chain files
REFERENCE_HG19=/genomes/scratch/dkasperaviciute/conversion/input/hg19.fa
REFERENCE_HG38=/genomes/scratch/dkasperaviciute/conversion/input/hg38.fa
REFERENCE_G37=/genomes/scratch/dkasperaviciute/conversion/input/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

CHAIN_G37_HG19=/genomes/software/apps/LiftOver/GRCh37ToHg19.over.chain.gz
CHAIN_HG19_HG38=/genomes/software/apps/LiftOver/hg19ToHg38.over.chain.gz
CHAIN_HG38_HG19=/genomes/software/apps/LiftOver/hg38ToHg19.over.chain.gz
CHAIN_HG19_G37=/genomes/software/apps/LiftOver/hg19ToGRCh37.over.chain.gz

INPUT_VCF=/home/mzarowiecki/scratch/TRACERx/3_external_exomes/Chosen2/LP3000396-DNA_H02.vcf.sort.vcf.gz
RESULTS_DIR=/home/mzarowiecki/scratch/TRACERx/3_external_exomes/2_liftOver2GRCh38/results/TRACERx/2018-06-19
PERL_HELPER=/home/mzarowiecki/git/AssemblyConversion/remove_diff_backconverted.pl
NORMALISE=/home/mzarowiecki/git/AssemblyConversion/normalise_vcf_dedup.sh

#INPUT_VCF=/genomes/resources/genomeref/data/human/GRCh37_annotation/gnomad-public/release-170228/vcf/genomes/gnomad.genomes.r2.0.1.sites.22.vcf.gz
#RESULTS_DIR=/genomes/scratch/dkasperaviciute/conversion/annotation/test2
#PERL_HELPER=/home/dkasperaviciute/git/AssemblyConversion/remove_diff_backconverted.pl

#mkdir -p $RESULTS_DIR
PREFIX=`basename $INPUT_VCF .vcf.gz`

#check for duplicated variants in input file
INPUT_DUPLICATES=`bcftools view -H $INPUT_VCF | cut -f 1,2,4,5 | sort | uniq -D | wc -l`
echo "`$NOW`INFO: number of duplicated records in input $PREFIX.vcf.gz $INPUT_DUPLICATES"

#GRCH37>hg19 conversion is required in order to work with UCSC chain files
#it should be perfect for chr 1-22, X and Y
#don't write original position here, will write in the next step, hg19 > hg38
echo "`$NOW`converting $PREFIX GRCh37 > hg19..."

picard LiftoverVcf \
	I=$INPUT_VCF \
	O=$RESULTS_DIR/$PREFIX.picard.hg19.vcf \
	CHAIN=$CHAIN_G37_HG19 \
	REJECT=$RESULTS_DIR/$PREFIX.picard.hg19.rejected.vcf \
	R=$REFERENCE_HG19

echo "`$NOW`converting $PREFIX hg19 > hg38..."

picard LiftoverVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg19.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg38.unsorted.vcf \
	CHAIN=$CHAIN_HG19_HG38 \
	REJECT=$RESULTS_DIR/$PREFIX.picard.hg38.rejected.vcf \
	R=$REFERENCE_HG38 \
	WRITE_ORIGINAL_POSITION=true

#sorting just in case

picard SortVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg38.unsorted.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg38.vcf

rm $RESULTS_DIR/$PREFIX.picard.hg38.unsorted.vcf

echo "`$NOW`removing entries with identical positions and alleles..."
#need to remove both occurencies of duplicated entries
#find line numbers of duplicated entries
#picard outputs a sorted file, therefore uniq should work without additional sorting

nl $RESULTS_DIR/$PREFIX.picard.hg38.vcf | cut -f 1,2,3,5,6 |  uniq -f 1 -D | cut -f 1 > $RESULTS_DIR/$PREFIX.duplicated.lines.hg38.txt

REMOVED_DUP=`wc -l $RESULTS_DIR/$PREFIX.duplicated.lines.hg38.txt`
echo "`$NOW`INFO: removed records converted to identical positions (number of lines removed): $REMOVED_DUP"

#now remove these lines from vcf only if file is not empty
if [[ -s $RESULTS_DIR/$PREFIX.duplicated.lines.hg38.txt ]] ; then
	sed 's%$%d%' $RESULTS_DIR/$PREFIX.duplicated.lines.hg38.txt | sed -f - $RESULTS_DIR/$PREFIX.picard.hg38.vcf > $RESULTS_DIR/$PREFIX.picard.hg38.no_dup.vcf
	mv $RESULTS_DIR/$PREFIX.picard.hg38.vcf $RESULTS_DIR/$PREFIX.picard.hg38.with_duplicates.vcf
	mv $RESULTS_DIR/$PREFIX.picard.hg38.no_dup.vcf $RESULTS_DIR/$PREFIX.picard.hg38.vcf
fi

echo "`$NOW`converting $PREFIX hg38 > hg19..."
#don't write original position, because need to preserve that field from earlier file

picard LiftoverVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg38.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg19.back.vcf \
	CHAIN=$CHAIN_HG38_HG19 \
	REJECT=$RESULTS_DIR/$PREFIX.picard.hg19.back.rejected.vcf \
	R=$REFERENCE_HG19

echo "`$NOW`removing records for variants converted to a different position from the original position"
#helper perl script
perl $PERL_HELPER $RESULTS_DIR/$PREFIX.picard.hg19.back.vcf $RESULTS_DIR/$PREFIX.picard.hg19.back.good.vcf $RESULTS_DIR/$PREFIX.picard.hg19.back.bad.vcf

echo "`$NOW`final conversion $PREFIX hg19 > hg38..."

picard LiftoverVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg19.back.good.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg38.final.unsorted.vcf \
	CHAIN=$CHAIN_HG19_HG38 \
	REJECT=$RESULTS_DIR/$PREFIX.picard.hg38.final.rejected.vcf \
	R=$REFERENCE_HG38

#sorting final file, just in case

picard SortVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg38.final.unsorted.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg38.final.vcf

rm $RESULTS_DIR/$PREFIX.picard.hg38.final.unsorted.vcf $RESULTS_DIR/$PREFIX.picard.hg38.final.unsorted.vcf.idx

#sanity check for duplicated positions
FINAL_DUPLICATES=`cut -f 1,2,4,5 $RESULTS_DIR/$PREFIX.picard.hg38.final.vcf | sort | uniq -D | wc -l`
echo "`$NOW`INFO: duplicated variants in final file $PREFIX.picard.hg38.final.vcf $FINAL_DUPLICATES"

#sanity check for back conversion

echo "`$NOW`sanity check back converting $PREFIX hg38 > hg19..."
#don't write original position, because need to preserve that field from earlier file

picard LiftoverVcf \
	I=$RESULTS_DIR/$PREFIX.picard.hg38.final.vcf \
	O=$RESULTS_DIR/$PREFIX.picard.hg19.back.sanity.vcf \
	CHAIN=$CHAIN_HG38_HG19 \
	REJECT=$RESULTS_DIR/$PREFIX.picard.hg19.back.sanity.rejected.vcf \
	R=$REFERENCE_HG19

echo "`$NOW`sanity check: removing records for variants converted to a different position from the original position"
#helper perl script
perl $PERL_HELPER $RESULTS_DIR/$PREFIX.picard.hg19.back.sanity.vcf $RESULTS_DIR/$PREFIX.picard.hg19.back.good.sanity.vcf $RESULTS_DIR/$PREFIX.picard.hg19.back.bad.sanity.vcf

echo "`$NOW`bgzip and index all vcf files"

for VCF in `find $RESULTS_DIR -name "$PREFIX.*.vcf"`; do
	echo "`$NOW`bgzip and index $VCF..."
	bgzip $VCF
	bcftools index $VCF.gz
done

echo "`$NOW`running normalistion"

$NORMALISE $RESULTS_DIR $PREFIX.picard.hg38.final $RESULTS_DIR/$PREFIX.picard.hg38.final.vcf.gz $REFERENCE_HG38

echo "`$NOW`done"



