#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAC-seq mapping pipeline (written by Charlie George)
Exercise 3/2/21 working directory /ifs/obds-training/jan21/exercises/peakcalling_mapping/01_pipeline_run_full_data

input:
    - Paired end ATACseq fastq.gz files & thier associated input files

outputs:
    - bamfiles
    - multiqc html report

processing:
    - fastqc
    - trimming
    - map reads using bowtie2

### TODO
### - make compatible with SE reads
### - parameterise trimming step
"""
#import required tools
import sys
from ruffus import *
from cgatcore import pipeline as P

#Define paramaters dictionary
params = P.get_parameters("atacseq_pipeline.yml")

#make directory called fastqc
@follows(mkdir("fastqc"))
#do fastqc on paired end reads
@transform("*.fastq.*.gz",
            regex(r"(.*).fastq.(.*).gz"),
            r"fastqc/\1_fastqc.\2.html")
def fastqc(infile, outfile):
    ''' run fastqc on all fastq.gz files

    input:
        string representing gzipped fastq file name e.g. sampleA.fastq.1.gz

    output:
        string representing fastqc html report name for fastq file e.g. fastqc/sampleA.fastq.1.gz

    example_statmenet:
        fastqc --nogroup -o fastqc sampleA.fastq.1.gz > fastqc/sampleA.fastq.1.gz.log
     '''

    statement = "fastqc --nogroup -o fastqc %(infile)s > %(outfile)s.log"
    P.run(statement, job_queue=params["queue"])


#make directory called reports
@follows(mkdir("reports"))
#run multiqc on fastqc output
@merge(fastqc, "reports/fastqc_report.html")
def multiqc(infiles, outfile):
    ''' run multiqc to collect fastq stats
    input:
        string of all files output from fastqc
    output:
        string of multiqc file html report
        reports/fastqc_report.html

    example statement:
        export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&
        multiqc fastqc/ -f -n reports/fastqc_report.html
    '''

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -f -n %(outfile)s"""
    P.run(statement, job_queue=params["queue"])


#make directory called trimmed_fastqs
@follows(mkdir("trimmed_fastqs"))
#perform trimming to remove adaptors from paired-end reads
@transform("*.fastq.1.gz",
            regex(r"(.*).fastq.1.gz"),
            r"trimmed_fastqs/\1.trimmed.fastq.1.gz")
def read_trimming(infile, outfile):
    ''' trim adapters from paired_end reads reads '''

    # get infile and outfile name for read1
    infile_read1 = infile
    outfile_read1 = outfile

    # create a sample name and outfile name without the read info on the end
    # e.g. sampleA.fastq.1.gz -> sampleA.fastq
    sample_name = P.snip(infile,'.1.gz')
    outfile_name = P.snip(outfile,'.1.gz')

    # create name of infile and outfile for read2 by adding '.2.gz' to end
    infile_read2 = f'{sample_name}.2.gz'
    outfile_read2 = f'{outfile_name}.2.gz'

    statement = '''trimmomatic PE
                        -threads 4
                        -phred33
                        %(infile_read1)s %(infile_read2)s
                        %(outfile_read1)s %(outfile_read1)s.unpaired
                        %(outfile_read2)s %(outfile_read2)s.unpaired
                        ILLUMINACLIP:/ifs/apps/bio/trimmomatic-0.32/adapters/NexteraPE-PE.fa:1:30:10:1:True
                        LEADING:3
                        TRAILING:3
                        SLIDINGWINDOW:4:15
                        MINLEN:30
                        2>> %(outfile)s.trimmomatic.log'''
    P.run(statement, job_queue=params["queue"],job_threads=4)



@merge(read_trimming, "reports/fastqc_trimming_report.html")
def multiqc_trimming(infiles, outfile):
    ''' run multiqc to collect fastq trimming stats '''

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc . -f -n %(outfile)s"""
    P.run(statement, job_queue=params["queue"])


@follows(multiqc_trimming,
         read_trimming)
@follows(mkdir("bam"))
@collate("trimmed_fastqs/*.trimmed.fastq.*.gz",
        regex(r"trimmed_fastqs/(.+).fastq.[1-2].gz"),
        r"bam/\1.bam")
def bowtie2(infiles, outfile):
    '''align reads to genome using bowtie2'''

    # infiles = tuple of read1 & read2 names
    # need to unpack it
    infile_read1 = infiles[0]
    infile_read2 = infiles[1]

    statement = """bowtie2 %(bowtie2_options)s
                        -p %(bowtie2_threads)s
                        -x %(bowtie2_index)s
                        -1 %(infile_read1)s
                        -2 %(infile_read2)s
                         2> %(outfile)s.log |
                   samtools sort -o %(outfile)s -@ %(bowtie2_threads)s -
                   && samtools index %(outfile)s"""
    P.run(statement, job_queue=params["queue"],
                     job_threads=params["bowtie2_threads"],
                     job_memory=params["bowtie2_memory"])

@transform(bowtie2,
            regex(r'(.*).bam'),
            r'\1.picardmetrics')
def picard(infile, outfile):
    ''' run picard to get alignment summary metrics'''
#set final memory for cluster to be slightly greater than memory used by picard tools
    final_memory = str(int(params['picard_memory']) + 2) + 'g'

    # VALIDATION_STRINGENCY=SILENT because getting error that it didn't like @PG header line?
    statement = """picard -Xmx%(picard_memory)sg
                            CollectAlignmentSummaryMetrics
                            R=%(picard_genome)s
                            I=%(infile)s O=%(outfile)s
                            VALIDATION_STRINGENCY=SILENT
                            """
    P.run(statement, job_queue=params['queue'], job_memory=final_memory)


@transform(bowtie2,
            regex(r'(.*).bam'),
            r'\1.idxstats')
def idxstats(infile, outfile):
    '''get idxstats'''

    statement = """samtools idxstats %(infile)s > %(outfile)s"""
    P.run(statement, job_queue=params['queue'])


@transform(bowtie2,
            regex(r'(.*).bam'),
            r'\1.flagstats')
def flagstats(infile, outfile):
    '''get flagstats'''

    statement = """samtools flagstat %(infile)s > %(outfile)s"""
    P.run(statement, job_queue=params['queue'])


@follows(picard,
         flagstats,
         idxstats)
@merge(flagstats, "reports/bam_report.html")
def multiqc3(infiles, outfile):
    ''' run multiqc to collect fastq stats '''

    statement = """export LC_ALL=en_US.UTF-8;
                   export LANG=en_US.UTF-8;
                   multiqc . -f -n %(outfile)s"""
    P.run(statement, job_queue=params["queue"])


@mkdir('test_bams')
@transform(bowtie2,
            regex(r'(.*).bam'),
            r'test_bams\1.chr5.bam')
def subset_bam(infile, outfile):
    ''' extract smaller version of bam '''

    statement = """samtools view -b %(infile)s chr5 > %(outfile)s.unsorted &&
                   samtools sort -o %(outfile)s -@ %(bowtie2_threads)s %(outfile)s.unsorted
                   && samtools index %(outfile)s"""
    P.run(statement, job_queue=params['queue'], job_threads=params["bowtie2_threads"])


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
