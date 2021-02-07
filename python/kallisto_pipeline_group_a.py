""" this is the python script for the alignment exercises on 2nd Feb 2021"""
""" Write a pseudoalignment pipeline"""
""" Use atom to edit files on the server"""
""" check your code into the shared github repo

1. Run fastqc on fastq.gz RNA-seq files (same as yesterday)
   Usage: fastqc -o output_directory --nogroup input_files
   What does nogroup mean? default for long reads
2. Run multiqc on fastqc output to generate a single report (same as yesterday)
   Useage: multiqc -n output_name.html -o output_directory input_directory
3. Index the mouse reference cDNA fasta using kallisto index
   Useage: kallisto index -i output_file input_file
4. Quantify the numbe of RNA-seq reads in the fastq.gz files using kallisto quant
  (refering to the indexed mouse reference). Save the stdout in a log file to be
   used in the next step.
   Usage: kallisto quant -i index [mouse file] -o output_directory [kallisto_quant]
   -b bootstraps -t threads pairA_1.fastq [input file 1] pairA_2.fastq [input file 2]
5. Run multiqc on the kallisto log file to generate a report in html
   Usage: multiqc -m kallisto -n output_file.html -o output_directory input_directory"""

#import required packages
from ruffus import *
from cgatcore import pipeline as P
import sys
import os

#Set parameters using the accompanying .yml file
PARAMS = P.get_parameters("kallisto_pipeline_group_a.yml")

#Make a new directory called fastqc
@follows (mkdir ("fastqc"))
#Run fastqc on the fastq.gz files
#the regex command switches the name of the input file from file.fasta.gz to file_fastqc.html
#> %(outfile)s.log 2>&1 sends stdout and stderr to log file
@transform ("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html')
def fastqc (infile, outfile):
    statement = '''fastqc -o fastqc --nogroup %(infile)s > %(outfile)s.log 2>&1'''
    P.run(statement)

#Run multiqc on the fastqc output
#%(outfile)s.log 2>&1 directs the stdout and stderr to a log file
@merge (fastqc, "fastqc/multiqc_fastqc.html")
def multiqc_fastqc (infiles, outfile):
    statement = '''multiqc -n multiqc_fastqc.html -o fastqc fastqc > %(outfile)s.log 2>&1'''
    P.run(statement)

#Make a directory called index (in retrospect this made later steps more complex)
@follows (mkdir ("index"))
#index the fasta file of mouse cDNA, which is listed in the .yml file
@transform(PARAMS["kallisto_index_genome_fasta"], regex(r'(.*).fa.gz'), r'index/\1.idx')
def kallisto_index (infile, outfile):
    statement = """kallisto index -i %(outfile)s  %(infile)s"""
    P.run(statement)


#Make a directory called kallisto_quant
@follows (mkdir ("kallisto_quant"))
#Quantify the number of reads using kallisto quant
@collate ("*.fastq.gz", regex(r'(.*)_[12].fastq.gz'), r'kallisto_quant/\1/abundance.tsv')
def kallisto_quant (infiles, outfile):
    #Pair the two input files for each sample (_1 and _2)
    infile1, infile2 = infiles
    #Specify the indeex
    index = os.path.basename(PARAMS["kallisto_index_genome_fasta"]).replace(".fa.gz", ".idx")
    #Make output folders for the results according to the input file name
    #(since output file names are fixed in kallisto)
    output_folder = outfile.replace("abundances.tsv", "")
    #Run kallisto quant with 10 bootstraps (from .yml) and multi threads
    #Send stdout and stderr to a log file
    statement = """kallisto quant -i index/%(index)s -o %(output_folder)s
                    -t %(kallisto_quant_threads)s -b %(kallisto_quant_bootstraps)s %(kallisto_quant_options)s
                    %(infile1)s %(infile2)s > %(outfile)s.log 2>&1"""
    #Thread info also needs to go in the P.run line
    P.run(statement, job_threads = PARAMS["kallisto_quant_threads"])

# run the multiqc on the output of kallisto_quant (abundances.tsv.log output)
@merge (kallisto_quant, "kallisto_quant/multiqc_kallisto.html")
def multiqc_abundance (infile, outfile):
    statement = '''multiqc -m kallisto -n multiqc_kallisto.html -o kallisto_quant
    kallisto_quant > %(outfile)s.log 2>&1'''
    P.run(statement)

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
