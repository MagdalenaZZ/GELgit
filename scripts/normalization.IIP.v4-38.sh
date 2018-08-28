#!/bin/bash

# request Bourne shell as shell for job
#$ -S /bin/bash

# request the all.q queue
#$ -q all.q

# error and output go to the same location
#$ -j y

#$ -cwd

module load bcftools/1.3
module load htslib/1.3
module load vt

# reference sequence is needed to left-align indels
if [ $REF_VERSION == "GRCh37" ]; then
	REFERENCE=/genomes/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa
elif [ $REF_VERSION == "GRCh38Decoy" ]; then
	REFERENCE=/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa
fi

VCF=$(awk "NR==$SGE_TASK_ID" $VCF_DIR/list/tumour_genomes.$TODAY.txt)
VCF_PREF=`echo "${VCF%%.*}"`
VCF_PASS_ONLY=${VCF_PREF}.PASS.somatic.vcf.gz
echo "input file ${VCF_PASS_ONLY}"

FILTERS=`gzip -c -d $VCF|grep FILTER|cut -s -d '=' -f 3|cut -d ',' -f 1|grep -v PASS`
FILTERS=`echo "$FILTERS"|perl -e '$filters = ""; while(<>){chomp(); $filters.="$_"."\|"} chop($filters) ;print "$filters";'`
gzip -c -d $VCF|grep -vP "$FILTERS"|gzip -c - > $VCF_PASS_ONLY

# split records with multiple alternate alleles into multiple biallelic records
# an e.g. 1/2 genotype will be split to 1/. and ./1
# You must specify the -s ("smart") option in order for INFO and FORMAT fields of type A and R to be retained and decomposed appropriately.
echo "###split records with multiple alternate alleles into multiple biallelic records"
vt decompose -s -o ${VCF_PREF}.PASS.split.somatic.vcf.gz $VCF_PASS_ONLY
bcftools index ${VCF_PREF}.PASS.split.somatic.vcf.gz 

# left-align indels and trim redundant bases
# check for consistency with the reference sequence
echo "###left-align indels and trim redundant bases"
vt normalize -r $REFERENCE -o ${VCF_PREF}.PASS.left.split.somatic.vcf.gz ${VCF_PREF}.PASS.split.somatic.vcf.gz 
bcftools index ${VCF_PREF}.PASS.left.split.somatic.vcf.gz 

# decompose complex variants to allelic primitives
# this will only operate on variants where REF and ALT alleles have the same length unless the -a option is specified
# not specifying -a for now
echo "###decompose complex variants to allelic primitives"
vt decompose_blocksub -o ${VCF_PREF}.PASS.atomic.left.split.somatic.vcf.gz ${VCF_PREF}.PASS.left.split.somatic.vcf.gz
bcftools index ${VCF_PREF}.PASS.atomic.left.split.somatic.vcf.gz 

# remove duplicates
# this program requires all fields of records to be identical before they are considered duplicates
# cf vt uniq, which defines duplicates based only on CHROM, POS, REF, ALT
echo "###remove duplicates"
bcftools norm -d both -Oz -o ${VCF_PREF}.PASS.duprem.atomic.left.split.somatic.vcf.gz ${VCF_PREF}.PASS.atomic.left.split.somatic.vcf.gz
tabix ${VCF_PREF}.PASS.duprem.atomic.left.split.somatic.vcf.gz

echo "###run GATK ValidateVariants"
java -jar /genomes/software/apps/gatk/3.4-46/GenomeAnalysisTK.jar \
   -T ValidateVariants \
   -R $REFERENCE \
   -V ${VCF_PREF}.PASS.duprem.atomic.left.split.somatic.vcf.gz \
   --validationTypeToExclude ALL
   

rm ${VCF_PREF}.PASS.somatic.vcf.gz*   
rm ${VCF_PREF}.PASS.split.somatic.vcf.gz*
rm ${VCF_PREF}.PASS.left.split.somatic.vcf.gz*
rm ${VCF_PREF}.PASS.atomic.left.split.somatic.vcf.gz* 

