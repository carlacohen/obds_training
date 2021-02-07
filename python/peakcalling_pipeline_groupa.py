""" groupa peakcalling_pipeline_groupa
Make a peakcalling pipeline with 4 tasks
Start from ATAC-seq reads which were trimmed and aligned to genome in previous pipeline.
So we have one file per sample now.
Using subset of test files comprising chr5 only (for quicker analysis)
4/2/21 working directory /ifs/obds-training/jan21/exercises/peakcalling/groupa
Contains softlinks to bam data files and this python script and associated .yml

1. Remove duplicate reads using picard mark_duplicates, also index the bams
Usage: picard MarkDuplicates #which bit of picard you want
                -Xmx%(picard_memory)sg # set memory use
                REMOVE_DUPLICATES=true # we are removing duplicates
                I=input_file (bam file from previous pipline i.e. RNAseq data that has been mapped to genome)
                O=output_file (bam file with duplicates removed)
                M=%(metrics)s.metrics
                VALIDATION_STRINGENCY=SILENT
                > %(outfile)s.log 2>&1 &&
                samtools index %(outfile)s #index the bam files with samtools
Set the queue and memory use in the P.run line.
2. Remove unmapped, unpaired, multimapping alignments and reads on chrM using samtools
Usage: samtools view %(samtools_flags)s #select flags corresponding to unmapped (0x4) and unpaired (0x2) reads
                -q 40 # select MAPQ score of >40 (identify uniquely mapped reads)
                -b input_file
                [select regions] #to exclude chrM mitochondrial reads
                > output
                && samtools index %(outfile)s #index the bam file
3. Remove reads from bam file in blacklist regions using bedtools intersect
Usage: bedtools intersect -v #-v means exclude selected regions
                -a input #bam from previous step
                -b input blacklist file from yml
                > output_file
                && samtools index output_file #also do index using samtools
4. Call peaks using Macs2
Set new conda environment using job_condaenv = "peaktools_env" in the P.run line.
Usage: macs2 callpeak
    --format BAMPE #This means Bam file from paired end reads
    --treatment input_file
    --verbose 10 # how much info to receive
    --bdg --SPMR --call-summits --keep-dup all #various other setting we had to look up (see below)
                                            #these could be parameterised
    --name output_file
    --gsize mm #mouse genome size
    >& %(outfile)s.log #send stdout to log file

Fasta file for picard (reference genome)
/ifs/mirror/genomes/plain/mm10.fasta ?do we ever use this?

Blacklist (reference)
/ifs/obds-training/jan21/exercises/peakcalling/data/proatac_combined_greenleaf_encode_blacklist_mm10.bed

"""
#Import required tools
import sys
from ruffus import *
from cgatcore import pipeline as P

#Set file containing parameters
PARAMS = P.get_parameters("peakcalling_pipeline_groupa.yml")

# make a directory called dedup_bams
@follows(mkdir("dedup_bams"))
# use a transform to convert one bam file into one deduplicated bam file
@transform("*.bam",
                regex(r"(.*).bam"),
                r"dedup_bams/\1.dedup.bam")
def picard_dedup (infile, outfile):
    metrics = outfile.replace(".bam", ".metrics")
    final_memory = str(int(PARAMS['picard_memory']) + 2) + 'g'
    statement = """picard MarkDuplicates
                    -Xmx%(picard_memory)sg
                    REMOVE_DUPLICATES=true
                    I=%(infile)s
                    O=%(outfile)s
                    M=%(metrics)s.metrics
                    VALIDATION_STRINGENCY=SILENT
                    > %(outfile)s.log 2>&1 &&
                    samtools index %(outfile)s """
    P.run(statement, job_queue=PARAMS['queue'], job_memory = final_memory)

#make a directory called samtools_filter
@follows(mkdir("samtools_filter"))
#Use samtools to filter the deduplicated bam files
@transform(picard_dedup,
                regex(r"dedup_bams/(.*).dedup.bam"),
                r"samtools_filter/\1.filter.bam")
def samtools_filter (infile, outfile):
    statement = """samtools view
                    %(samtools_flags)s
                    -q 40
                    -b %(infile)s
                    %(samtools_chromosomes)s
                    > %(outfile)s
                    && samtools index %(outfile)s"""
    P.run(statement)

#make directory called bed_blist_intersect
@follows(mkdir('bed_blist_intersect'))
# use bedtools intersect to find NON-OVERLAPPING regions with blacklist
@transform(samtools_filter,
            regex(r"samtools_filter/(.*).filter.bam"),
            r"bed_blist_intersect/\1.filter.intersect.bam")
def bed_blist_intersect (infile, outfile):
    statement = """bedtools intersect -v -a %(infile)s -b %(bedtools_blacklist)s
        > %(outfile)s
        && samtools index %(outfile)s"""
    P.run(statement)

#make directory called macs2_peak
@follows(mkdir('macs2_peak'))
#run macs2 callpeak
@transform(bed_blist_intersect,
            regex(r"bed_blist_intersect/(.*).intersect.bam"),
            r"macs2_peak/\1.macs2_peak.bam")
def macs2_peak (infile, outfile):
    statement = """macs2 callpeak
        --format BAMPE
        --treatment %(infile)s
        --verbose 10
        --bdg --SPMR --call-summits --keep-dup all
        --name %(outfile)s
        --gsize mm
        >& %(outfile)s.log"""
    # new conda env needed that contains macs2
    P.run(statement, job_condaenv = "peaktools_env")

# --bdg is the number of fragments within the pileup
# --SPMR needed to visualise in genome browser
# --call-summits is recommended to detect adjucent binding events/distinguishes two nearby peaks
# --tsize means size of sequencing tags

#The if__name__ == "__main__": Whenever the Python interpreter reads a source file, it does two things:
#it sets a few special variables like __name__, and then
#it executes all of the code found in the file.

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
