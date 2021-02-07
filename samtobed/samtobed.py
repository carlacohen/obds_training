#! /usr/bin/env python
# allow script to run on command line (doesn't work)
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:46:36 2021

Aims of the exercise:
1. Write a Python script to convert the SAM file to a BED file
• One BED line per SAM line
• Read SAM file input
• Write BED file output - format your output using fstrings

NB BED format is: chr start end score name
How do we get these from the sam file?
chr - comes from 3rd column 
start - is column 4
end - need to calculate the length of each sequence from column 10, then can calculate end
name - is column 1
score - use mapping quality score from column 5
strand - use + for all (mapping direction)

Need to remove the header from the file which starts with @

2. Supply input file names using command line parameters
• Supply the SAM file name on the command line using –i or --input
• Supply the BED file name on the command line using –o or --output

3. Write out a gzip compressed file (.bed.gz)

4. Make the script run on the command line

@author: ccohen
"""

#import the packages that you need
import argparse
import gzip

#Set the parameters
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", dest="bedfilename", default="output.bed.gz", help="output file name")
parser.add_argument("-i", "--input", dest="samfilename", help="Input file name")
parser.add_argument("-p", "--pad", dest="padvalue", default=0, help="Input pad number to extend start and end position")
parser.add_argument("-f", "--fragment", dest="fragment", default=False, action='store_true', help="Output frgament coordinates rather than read coordinates")
args = parser.parse_args()

#create a bed.gz file to write into
bedfile = gzip.open(args.bedfilename, "wt")  

# Open a file for reading
with open (args.samfilename, "r") as samfile:
    #iterate line by line
    for line in samfile:
        #remove header lines starting with @
        if not line.startswith("@"):
            #separate the columns by tab
            columns = line.split("\t")
            #select the right column for each bed column
            chrom = columns[2]
            #defining insert size
            insert_size = int(columns[8])
            #ignore unmapped reads with * in column 3
            if not chrom == "*":    
                start = int(columns[3]) - int(args.padvalue)
            #calculate the end and tell it that col3 is an integer
                if args.fragment: 
                    end = start+insert_size + int(args.padvalue)
                else:
                    end = len(columns[9])+int(columns[3]) + int(args.padvalue)
                name = columns[0]
                score = columns[4]
                #write into the bedfile
                if args.fragment: 
                    if insert_size >0: 
                        bedfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t+\n")
                else:
                    bedfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t+\n")        

bedfile.close()
  