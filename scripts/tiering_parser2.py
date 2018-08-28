#!/usr/bin/env python
from __future__ import print_function
<<<<<<< HEAD

#import pysam
=======
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
import sys
import os.path
import argparse
import json

"""

Script for parsing and editing a tiering json file

"""


epi = ('\
    \n\
	json file parser, allowing for custom advanced filtering of tiering files\n\
    \n\
')


# Describe what the script does
<<<<<<< HEAD
parser = argparse.ArgumentParser(description='This script parses a json file and applies custom filtering', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='jso', action='store', required=True, help="Json file")
#parser.add_argument('-f', '--filter-string', default=None, dest='fi', action='store', required=True, help="String of filters to apply")
=======
parser = argparse.ArgumentParser(description='This script parses a json tiering file and applies custom filtering', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='jso', action='store', required=True, help="Json file")
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.jso)==True:
<<<<<<< HEAD
    print("Cannot find input file ",args.jso)
=======
    print("Cannot find input file ", args.jso)
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
    sys.exit(1)

# Read input
inp = open (args.jso)
json_array = json.load(inp)

# Create outfile
output=args.jso+".ann"
<<<<<<< HEAD
print("Input: ",args.jso, "Output: " , output)
out = open(output, "w")


#store_list = []
=======
output2=args.jso+".vars"
print("Input: ",args.jso, "Output: " , output)
out = open(output, "w")
out2 = open(output2, "w")

results = {}
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1

for j in json_array:
    chr = "chr" + j["reportedVariantCancer"]["chromosome"]
    pos = j["reportedVariantCancer"]["position"]
<<<<<<< HEAD

    dict = j["reportedVariantCancer"]["additionalTextualVariantAnnotations"]


    # Flag H Indels intersecting with reference homopolymers of at least 8 nucleotides in length
    if 'IHP' in dict:
        print(chr,pos,"H-Homopolymers", dict['IHP'], sep="\t", file = out)

    # Flag N	Small indels in regions with high levels of sequencing noise where at least 10% of the basecalls in a window extending 50 bases to either side of the indels call have been filtered out due to the poor quality
    if 'indel_noise' in dict:
        print(chr, pos, "N-SmallIndelsInNoise", dict['indel_noise'], sep="\t", file = out)
=======
    dict = j["reportedVariantCancer"]["additionalTextualVariantAnnotations"]
    vaf =  j["reportedVariantCancer"]["vaf"]
    pc = j["reportedVariantCancer"]["proteinChange"]



    key = '_'.join((chr,str(pos)))
    print (chr,str(pos), vaf, sep="\t", file = out2)

    # Check if some annotation has already been done, or add a new empty key to the results dict
    if key in results:
        pass
    else:
        results[key]={}

    # Flag H Indels intersecting with reference homopolymers of at least 8 nucleotides in length
    if 'IHP' in dict:
        #print(chr, pos, "H-Homopolymers", dict['IHP'], sep="\t", file = out)
        results[key]['H-Homopolymers']=dict['IHP']


    # Flag N	Small indels in regions with high levels of sequencing noise where at least 10% of the basecalls 
    # in a window extending 50 bases to either side of the indels call have been filtered out due to the poor quality
    if 'indel_noise' in dict:
        #print(chr, pos, "N-SmallIndelsInNoise", dict['indel_noise'], sep="\t", file = out)
        results[key]['N-SmallIndelsInNoise']=dict['indel_noise']

>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1

    # Flag G Variants with germline allele frequency > 1% in an internal Genomic England data set
    # (indicates potential un-subtracted germline variant)
    if 'GEL.GL.6628_AF' in dict:
        #print(chr, pos, "G-potential_germline", dict['GEL.GL.6628_AF'],"obsolete", sep="\t")
        # And the value is G
        if dict['GEL.GL.6628_AF']=='G':
<<<<<<< HEAD
            print(chr, pos, "G-potential_germline", dict['GEL.GL.6628_AF'], sep="\t", file = out)
=======
            #print(chr, pos, "G-potential_germline", dict['GEL.GL.6628_AF'], sep="\t", file = out)
            results[key]['G-potential_germline']=dict['GEL.GL.6628_AF']
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1


    # Flag R Recurrently identified somatic variants with somatic allele frequency > 5% in an
    # internal Genomic England data set (indicates potential technical artefact)
    #"somatic_agg_vcf_AF_FFnano": "R",
    #"somatic_agg_vcf_AF_FFpcrfree": "R",
    #"somatic_agg_vcf_AF_FFPE": "R",
<<<<<<< HEAD

    if 'somatic_agg_vcf_AF_FFnano' in dict:
        print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFnano'], sep="\t", file = out)
    if 'somatic_agg_vcf_AF_FFpcrfree' in dict:
        print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFpcrfree'], sep="\t", file = out)
    if 'somatic_agg_vcf_AF_FFPE' in dict:
        print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFPE'], sep="\t", file = out)
=======
    if 'somatic_agg_vcf_AF_FFnano' in dict:
        #print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFnano'], sep="\t", file = out)
        results[key]['R-RecurrentSomatic']=dict['somatic_agg_vcf_AF_FFnano']
    if 'somatic_agg_vcf_AF_FFpcrfree' in dict:
        #print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFpcrfree'], sep="\t", file = out)
        results[key]['R-RecurrentSomatic']=dict['somatic_agg_vcf_AF_FFpcrfree']
    if 'somatic_agg_vcf_AF_FFPE' in dict:
        #print(chr, pos, "R-RecurrentSomatic", dict['somatic_agg_vcf_AF_FFPE'], sep="\t", file = out)
        results[key]['R-RecurrentSomatic']=dict['somatic_agg_vcf_AF_FFPE']
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1


    # Flag SR Variants overlapping simple repeats
    if 'simple_repeat' in dict:
<<<<<<< HEAD
        print(chr, pos, "SR-SimpleRepeat", dict['simple_repeat'], sep="\t", file = out)

=======
        #print(chr, pos, "SR-SimpleRepeat", dict['simple_repeat'], sep="\t", file = out)
        results[key]['SR-SimpleRepeat']=dict['simple_repeat']


    
    #if 'vaf' in dict2:
    #    print(chr, pos, "vaf:", vaf, "pc:", pc , sep="\t", file = out2)



#for key in sorted(mydict.iterkeys()):
for chp, flag in results.items():

    if bool(flag)==True:
        print(chp.split("_")[0],"\t",chp.split("_")[1],"\t" , sep="", end="", file = out)
        print(','.join(flag), file = out)
        
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1

out.close()


<<<<<<< HEAD
quit()
=======
quit()
>>>>>>> 33eaa41eef46817145ed7ea3679e0c6b4a6271c1
