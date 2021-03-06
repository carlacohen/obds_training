'''Aim to write a bedtools pipeline to intersect ATAC and ChIP masterlists
with Refseq gene TSS
Author: Carla Cohen 7.2.21
Working directory /well/jknight/Carla/Python_test/bedtools

Input files:
    ATAC & ChIP masterlists, already sorted bedfiles
    Refseq TSS sorted bedfile

Output files:
    Text files of regions in ML bed that intersect within specified distance of TSS
    NB by adding %(distance)s.txt to the outfile name you get the distance in the file name for easy reference.

Usage: bedtools window -w [distance] -a input_ML -b refseq_genes_TSS.bed'''

# import packages
import sys
from ruffus import *
from cgatcore import pipeline as P

# read in paramerters from yml file
params = P.get_parameters("bedtools_pipeline.yml")

@mkdir("outputs")
@transform('*_chr5.bed',
            regex(r'(.*)_chr5.bed'),
            r'outputs/\1_chr5_TSS')
def bedtools_window(infile, outfile):
    statement = """bedtools window -w %(distance)s
                -a %(infile)s
                -b %(refseq_TSS)s
                > %(outfile)s.%(distance)s.txt"""
    P.run(statement)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
