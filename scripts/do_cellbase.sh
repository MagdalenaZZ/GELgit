#!/bin/bash

module load cellbase/v4.5.5

cellbase.sh variant-annotation \
 --input-file input.vcf.gz \
 --local \
 --assembly GRch38 \
 --output output.json.gz \
 --species hsapiens \
 --skip-decompose \
 --exclude expression,geneDisease,drugInteraction,conservation,functionalScore \
  --custom-file /genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/index_somatic_variants/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz,/genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/annotate_germline_variants/GEL_GL_6628.duprem.sites.annot.subset.atomic.left.split.vcf.gz \
 --custom-file-id somatic_agg_vcf,GEL.GL.6628 \
  --custom-file-fields AF_FFPE,AF_FFpcrfree,AF_FFnano:GN,AF,AC,AN,MAF,HWE,AN_Cancer,AN_SRv3,AN_RD,AN_SRv4,AC_Cancer,AC_SRv3,AC_RD,AC_SRv4,AF_Cancer,AF_SRv3,AF_RD,AF_SRv4,MAF_Cancer,MAF_SRv3,MAF_RD,MAF_SRv4,HWE_Cancer,HWE_SRv3,HWE_RD,HWE_SRv4
