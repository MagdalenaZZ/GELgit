#!/bin/bash

#usage sh codon_analysis_script.sh 7 140453135 140453145 <(samtools view ../RNA-Seq_Data/aligned/PADIMAC-119/PADIMAC-119_unique.bam 7:140453135-140453145)

#N.B. this will not print out codons where there is an insertion (or intron), which can mean the beginning is the other side of an intron, and so start + length will not be greater than codon end
chrom=$1
codon_start=$2
codon_end=$3

while read line; do
        start=$(echo "$line" | cut -f 4)
        sequence=$(echo "$line" | cut -f 10)
#        sequence=$(echo "$line" | cut -f 11)
        length=${#sequence}
#       codon=$(echo "$sequence" | cut -c 1-$(($start-$codon_start)))
        if [ $start -lt $codon_start ] && [ $(($start+$length)) -gt $codon_end ]
        then
#               echo "hello"
#               echo $(($codon_start-$start))
#               echo $(($codon_end-$codon_start))
                codon=$(echo "$sequence" | cut -c $(($codon_start-$start))-$(($codon_start-$start+$codon_end-$codon_start)))
                echo $codon
        fi
        #echo $codon
# | cut -c ($start+3)-($codon_end-$codon_start) $sequence | awk '{codon[$0]++} END{for(x in codon){print x,codon[x]}}';

done < $4
