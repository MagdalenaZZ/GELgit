#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import pysam
import sys
import os.path
import argparse

"""

Script for parsing and editing a VCF file

"""


epi = ('\
    \n\
	VCF file parser, allowing for custom advanced filtering of VCF files\n\
    \n\
')


# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and applies custom filtering', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', action='store', required=True, help="VCF file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.vcf)==True:
    print("Cannot find input file ",args.vcf)
    sys.exit(1)

<<<<<<< HEAD
=======
# Check if MNVs is to be split
if (args.sm)=="True":
    print ("MNVs will be split")
elif (args.sm) == "False":
    pass
    #print("MNVs will NOT be split")
else:
    print("MNVs option shoud be \"True\" or \"False\", is :", args.sm )

>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
# Read input
output=args.vcf+".pysam.vcf"
outs=args.vcf+".stats"



# read the input file
myvcf = pysam.VariantFile(args.vcf, "r")

myvcf.header.info.add("VT", "1", "String", "Variant type")

# create an object of new vcf file and open in to write data.
vcf_out = pysam.VariantFile(output, 'w', header=myvcf.header)
output=args.vcf+".sstats"

out = open(output, 'w')

# tx
totalx = 0  # tx
SNPx=0
indelx=0
otherx=0

# Tx VAF>10%
total10x = 0  # tx_vaf10
SNP10x=0
indel10x=0
other10x=0

# Tx VAF>20%
total20x = 0  # tx_vaf20
SNP20x=0
indel20x=0
other20x=0

# Our VAF>10%
total10 = 0  # vaf10
SNP10=0
indel10=0
other10=0

# Shared
shared_total = 0  # shared
shared_SNP=0
shared_indel=0
shared_other=0

# Shared + tx VAF>10%
shared_total10x = 0  # shared_tx_vaf10
shared_SNP10x=0
shared_indel10x=0
shared_other10x=0

# Shared + tx VAF>20%
shared_total20x = 0  # shared_tx_vaf10
shared_SNP20x=0
shared_indel20x=0
shared_other20x=0

# Shared + our VAF>10%
shared_total10 = 0  # shared_vaf10
shared_SNP10=0
shared_indel10=0
shared_other10=0

#
total=0
snp=0
indel=0
other=0


for r in myvcf:

    #### FILTER OUT #####
    # Shared called total
    # Filter out sites which
    refb = r.ref
    altb = r.alts
    #r.info["VT"] = "None"

    # Filter variants which are not in TRACERx
    if 'FQ' in r.samples[0].keys():
        if (''.join((r.samples[0]['FQ']))=='.'):
            r.info['VT'] = "Not_in_TRACERx"
            #print (''.join((r.samples[0]['FQ'])))
            continue

    # Filter TRACERx 0/0 gt
    if 'GT' in r.samples[0].keys():
        if (r.samples[0]['GT'][1]==0):
            r.info['VT'] = "TRACERx_0_GT"
            continue

    #### MARK/MODIFY #####

    tx = 0  # is a variant in tx
    tx_vaf10 = 0 # TX vaf > 10%
    tx_vaf20 = 0  # TX vaf > 10%

    vaf10 = 0  # our VAF > 10%

    shared = 0  # Shared pass TX + us
    shared_tx_vaf10=0 # shared + TX vaf > 10%
    shared_vaf10= 0  # shared + our VAF > 10%
    shared_tx_vaf20= 0  # shared + our VAF > 10%

    # Mark shared
    if 'QU' in r.samples[0].keys():

        # If TRACERx has a variant here any VAF
        if (r.samples[0]['QU'][0] == "PASS" or r.samples[0]['QU'][0] == "FAIL"):
            tx=1
        # If TRACERx has a variant here and VAF >10%
        if (r.samples[0]['QU'][0] == "PASS"):
            tx_vaf10 = 1
        if (r.samples[0]['QU'][0] == "PASS" and  (float(r.samples[0]['FQ'][0]) > 0.2) ):
            tx_vaf20 = 1
        # If TRACERx has a variant here, and we also have a PASS
        if (r.samples[1]['QU'][0] == "PASS"):
            shared = 1
        # If TRACERx has a variant here and VAF >10%, and we also have a PASS
        if (r.samples[0]['QU'][0] == "PASS" and r.samples[1]['QU'][0] == "PASS"):
            shared_tx_vaf10 = 1
            # If shared variant, and TRACERx VAF >20%
            if (float(r.samples[0]['FQ'][0]) > 0.2):
                #print("More than 20",r.samples[0]['QU'][0] , r.samples[1]['QU'][0] , r.samples[0]['FQ'][0] )
                shared_tx_vaf20 = 1

    # VAF more than 10%
    if 'FQ' in r.samples[1].keys():
        if (r.samples[1]['FQ'][0]=="."):
            #print("Dot ", r.samples[1]['FQ'][0])
            shared=0
            pass
        # If us have VAF > 10%
        elif (float(r.samples[1]['FQ'][0]) > 0.1):
            #print ("More than 0.1 ",r.samples[1]['FQ'][0])
            vaf10=1
            # If there is also a pass
            if (shared==1):
                shared_vaf10 = 1
        else:
            #print ("Less than 0.1 ",r.samples[1]['FQ'][0])
            pass


    # Filter multiple alts - indicating that we had different opinions what the variant was
    if (len(altb)>1):
        shared=0

        if (args.sm == "False"):
            r.info['VT'] = "Multiple_alts"
            continue

        # If MNVs are to be split, proceed
        if (args.sm=="True"):

            # Make all genotypes to lists
            ref=list(refb)
            alt0 = list(altb[0])
            alt1 = list(altb[1])

            print("Multiple alts:",r.ref, r.alts[0], r.alts[1])
            print (len(ref), len(alt0), len(alt1))

            # If they are all of equal length
            if (len(alt0) !=len(alt1)):
                print ("Not same length")
                continue
            else:
                print ("Same length")
        mulalts=mulalts+1


    #### TALLY UP #####

    total = total + 1
    if (tx==1):
        totalx = totalx + 1
    if (tx_vaf10 == 1):
        total10x = total10x + 1
    if (tx_vaf20 == 1):
        total20x = total20x + 1
    if (vaf10==1):
        total10 = total10 + 1
    if (shared==1):
        shared_total = shared_total + 1
    if (shared_tx_vaf10==1):
        shared_total10x = shared_total10x + 1
    if (shared_vaf10==1):
        shared_total10 = shared_total10 + 1
    if (shared_tx_vaf20==1):
        shared_total20x = shared_total20x + 1


    # Shared called SNPs, indel, other
    if (len(refb)==1 and len(altb[0])==1):
        print(r.chrom, r.pos, refb, altb[0], shared, "SNV", ':'.join((str(r.samples[0]['GT'][0]), str(r.samples[0]['GT'][1]) , str(r.samples[0]['FQ'][0]),  r.samples[0]['QU'][0]) ) , ':'.join((str(r.samples[1]['GT'][0]), str(r.samples[1]['GT'][1]) , str(r.samples[1]['FQ'][0]),  r.samples[1]['QU'][0]) ) , sep='\t' )
        #print ("SNP ",refb,altb[0])
        r.info['VT'] = "Snv"
        if (tx == 1):
            SNPx = SNPx + 1
        if (tx_vaf10 == 1):
            SNP10x = SNP10x + 1
        if (tx_vaf20 == 1):
            SNP20x = SNP20x + 1
        if (vaf10 == 1):
            SNP10 = SNP10 + 1
        if (shared == 1):
            shared_SNP = shared_SNP + 1
        if (shared_tx_vaf10 == 1):
            shared_SNP10x = shared_SNP10x + 1
        if (shared_vaf10 == 1):
            shared_SNP10 = shared_SNP10 + 1
        if (shared_tx_vaf20 == 1):
            shared_SNP20x = shared_SNP20x + 1

    # Shared called MNV
    elif (len(refb)==2 and len(altb[0])==2):
        print(r.chrom, r.pos, refb, altb[0], shared, "MNV", ':'.join((str(r.samples[0]['GT'][0]), str(r.samples[0]['GT'][1]) , str(r.samples[0]['FQ'][0]),  r.samples[0]['QU'][0]) ) , ':'.join((str(r.samples[1]['GT'][0]), str(r.samples[1]['GT'][1]) , str(r.samples[1]['FQ'][0]),  r.samples[1]['QU'][0]) ) , sep='\t' )
        #print("MNV ",refb, altb[0])
        r.info['VT'] = "Mnv"
        if (tx==1):
            otherx = otherx + 1
        if (tx_vaf10 == 1):
            other10x = other10x + 1
        if (tx_vaf20 == 1):
            other20x = other20x + 1
        if (vaf10==1):
            other10 = other10 + 1
        if (shared==1):
            shared_other = shared_other + 1
        if (shared_tx_vaf10==1):
            shared_other10x = shared_other10x + 1
        if (shared_vaf10==1):
            shared_other10 = shared_other10 + 1
        if (shared_tx_vaf20 == 1):
            shared_other20x = shared_other20x + 1


    # Shared called deletion
    elif (len(refb) > len(altb[0])):
        print(r.chrom, r.pos, refb, altb[0], shared, "Deletion", ':'.join((str(r.samples[0]['GT'][0]), str(r.samples[0]['GT'][1]) , str(r.samples[0]['FQ'][0]),  r.samples[0]['QU'][0]) ) , ':'.join((str(r.samples[1]['GT'][0]), str(r.samples[1]['GT'][1]) , str(r.samples[1]['FQ'][0]),  r.samples[1]['QU'][0]) ) , sep='\t')
        r.info['VT'] = "Deletion"
        if (tx==1):
            indelx = indelx + 1
        if (tx_vaf10 == 1):
            indel10x = indel10x + 1
        if (tx_vaf20 == 1):
            indel20x = indel20x + 1
        if (vaf10==1):
            indel10 = indel10 + 1
        if (shared==1):
            shared_indel = shared_indel + 1
        if (shared_tx_vaf10==1):
            shared_indel10x = shared_indel10x + 1
        if (shared_vaf10==1):
            shared_indel10 = shared_indel10 + 1
        if (shared_tx_vaf20 == 1):
            shared_indel20x = shared_indel20x + 1


    # Shared called insertion
    elif (len(refb) < len(altb[0])):
<<<<<<< HEAD
        print(r.chrom, r.pos, r.ref, altb[0],shared, "Insertion")
        #':'.join(str(r.samples[1]['GT'][0]),str(r.samples[1]['DP'][0]),str(r.samples[1]['FQ'][0]),str(r.samples[1]['QU'][0])))
        r.info['VT'] = "Indel"
        if (tx == 1):
=======
        print(r.chrom, r.pos, refb, altb[0], shared, "Insertion", ':'.join((str(r.samples[0]['GT'][0]), str(r.samples[0]['GT'][1]) , str(r.samples[0]['FQ'][0]),  r.samples[0]['QU'][0]) ) , ':'.join((str(r.samples[1]['GT'][0]), str(r.samples[1]['GT'][1]) , str(r.samples[1]['FQ'][0]),  r.samples[1]['QU'][0]) )  , sep='\t')
        r.info['VT'] = "Insertion"
        if (tx==1):
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
            indelx = indelx + 1
        if (tx_vaf10 == 1):
            indel10x = indel10x + 1
        if (tx_vaf20 == 1):
            indel20x = indel20x + 1
        if (vaf10 == 1):
            indel10 = indel10 + 1
        if (shared == 1):
            shared_indel = shared_indel + 1
        if (shared_tx_vaf10 == 1):
            shared_indel10x = shared_indel10x + 1
        if (shared_vaf10 == 1):
            shared_indel10 = shared_indel10 + 1
        if (shared_tx_vaf20 == 1):
            shared_indel20x = shared_indel20x + 1

    else:
        shared_other = shared_other + 1
        print("Else ", refb, altb[0])


    vcf_out.write(r)


## WRITE OUTPUT ##

out = open(outs, "w")


out.write(''.join(("TRACERx total:\t", str( totalx),"\n")))
out.write(''.join(("TRACERx SNP:\t", str( SNPx ),"\n")))
out.write(''.join(("TRACERx indels:\t", str(indelx ),"\n")))

out.write(''.join(("Total Tx VAF10:\t", str( total10x),"\n")))
out.write(''.join(("SNP Tx VAF10:\t", str(SNP10x ),"\n")))
out.write(''.join(("Indels Tx VAF10:\t", str(indel10x ),"\n")))

<<<<<<< HEAD
=======
out.write(''.join(("TRACERx total:\t", str(totalx),"\n")))
out.write(''.join(("TRACERx SNP:\t", str(SNPx),"\n")))
out.write(''.join(("TRACERx indels:\t", str(indelx ),"\n")))

out.write(''.join(("Total Tx VAF10:\t", str( total10x),"\n")))
out.write(''.join(("SNP Tx VAF10:\t", str(SNP10x ),"\n")))
out.write(''.join(("Indels Tx VAF10:\t", str(indel10x ),"\n")))

>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
out.write(''.join(("Total Tx VAF20:\t", str( total20x),"\n")))
out.write(''.join(("SNP Tx VAF20:\t", str(SNP20x ),"\n")))
out.write(''.join(("Indels Tx VAF20:\t", str(indel20x ),"\n")))

out.write(''.join(("Shared total:\t", str( shared_total),"\n")))
out.write(''.join(("Shared SNP:\t", str(shared_SNP),"\n")))
out.write(''.join(("Shared indels: \t", str(shared_indel ),"\n")))

out.write(''.join(("Shared total Tx VAF10:\t", str( shared_total10x),"\n")))
out.write(''.join(("Shared SNP Tx VAF10:\t", str(shared_SNP10x ),"\n")))
out.write(''.join(("Shared indels Tx VAF10:\t", str(shared_indel10x ),"\n")))

out.write(''.join(("Shared total Tx VAF20:\t", str( shared_total20x),"\n")))
out.write(''.join(("Shared SNP Tx VAF20:\t", str(shared_SNP20x ),"\n")))
out.write(''.join(("Shared indels Tx VAF20:\t", str(shared_indel20x ),"\n")))

out.write(''.join(("Shared total us VAF10:\t", str( shared_total10),"\n")))
out.write(''.join(("Shared SNP us VAF10:\t", str(shared_SNP10 ),"\n")))
out.write(''.join(("Shared indels us VAF10:\t", str(shared_indel10 ),"\n")))

out.write(''.join(("Us Total VAF10:\t", str( total10),"\n")))
out.write(''.join(("Us SNP VAF10:\t", str(SNP10 ),"\n")))
out.write(''.join(("Us Indels VAF10:\t", str(indel10 ),"\n")))


out.close()
<<<<<<< HEAD
=======
vcf_out.close()
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1



exit(0)

