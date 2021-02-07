'''
RNAseq pipeline

01/02/2021
Authors: Ian, Carla, Jian, Bjorn, Ahmed
Working directory Feb 2nd 2021 /ifs/obds-training/jan21/exercises/rnaseq_python
Contains softlinks to data RNAseq fastq data files, this script and associated yml.

Commands from week 1:
#Quality control of fastq files
1. fastqc -o fastqc --nogroup -t 2 ERR1755082_test_1.fastq.gz ERR1755082_test_2.fastq.gz
#Building a report over all the files
2. multiqc -n multiqc_fastqc.html -o fastqc fastqc
#Map reads to the reference genome
3. hisat2 --threads 12 -x /ifs/mirror/genomes/hisat2/mm10 -1 ERR1755082_1.fastq.gz
-2 ERR1755082_2.fastq.gz --rna-strandness RF --summary-file sam/ERR1755082.txt -S sam/ERR1755082.sam
#Convert sam to bam and sort and index
4. samtools sort -@ 8 -o ERR1755082.bam ERR1755082.sam && samtools index ERR1755082.bam
#Generate mapping QC
5. samtools idxstats ERR1755082.bam > ERR1755082.idxstats
6. samtools flagstat ERR1755082.bam > ERR1755082.flagstat
7. picard CollectAlignmentSummaryMetrics -I=input.bam O=output.picard -R /ifs/mirror/genomes/plain/mm10.fasta
8. picard CollectInsertSizeMetrics -H ERR1755082_histogram.pdf -I ERR1755082.bam -O ERR1755082_picard_insertsizemetrics &2 > ERR1755082_picard.log &
#Use multiqc to visualise the mapping
9. multiqc -n multiqc_sam.html .
#Use feature counts to generate a count matrix
10. featureCounts -a /ifs/obds-training/jan21/shared/linux/rnaseq.genes.gtf.gz -o counts/countdata.tsv
'''
from ruffus import *
from cgatcore import pipeline as P
import sys

PARAMS = P.get_parameters("rnaseq_pipeline.yml")

@follows (mkdir ("fastqc"))
@transform ("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html')
def fastqc (infile, outfile):
    statement = '''fastqc -o fastqc --nogroup %(infile)s > %(outfile)s.log 2>&1'''
    P.run(statement)

@merge (fastqc, "fastqc/multiqc_fastqc.html")
def multiqc_fastqc (infiles, outfile):
    statement = '''multiqc -n multiqc_fastqc.html -o fastqc fastqc > %(outfile)s.log 2>&1'''
    P.run(statement)

@follows (mkdir ("bam"))
@collate ("*.fastq.gz", regex(r'(.*)_[12].fastq.gz'), r'bam/\1.bam')
def hisat2 (infiles, outfile):
    infile1, infile2 = infiles
    statement = '''hisat2 --threads %(hisat2_threads)s -x %(hisat2_genome)s -1 %(infile1)s
    -2 %(infile2)s %(hisat2_options)s --summary-file %(outfile)s.log |
    samtools sort -@ %(hisat2_threads)s -o %(outfile)s - && samtools index %(outfile)s'''
    P.run(statement, job_threads = PARAMS["hisat2_threads"])

@transform (hisat2, regex(r'bam/(.*).bam'), r'bam/\1.idxstats')
def idxstats (infile, outfile):
    statement = '''samtools idxstats %(infile)s > %(outfile)s'''
    P.run(statement)

@transform (hisat2, regex(r'bam/(.*).bam'), r'bam/\1.flagstat')
def flagstat (infile, outfile):
    statement = '''samtools flagstat %(infile)s > %(outfile)s'''
    P.run(statement)

@transform (hisat2, regex(r'bam/(.*).bam'), r'bam/\1.picard_metrics')
def picard_metrics (infile, outfile):
    statement = '''picard CollectAlignmentSummaryMetrics -I %(infile)s -O %(outfile)s -R %(picard_genome)s 2> %(outfile)s.log '''
    P.run(statement)

@transform (hisat2, regex(r'bam/(.*).bam'), r'bam/\1.picard_isize')
def picard_isize (infile, outfile):
    statement = '''picard CollectInsertSizeMetrics -I %(infile)s -O %(outfile)s -R %(picard_genome)s
    -H %(outfile)s.pdf 2> %(outfile)s.log '''
    P.run(statement)

@merge ([idxstats, flagstat, picard_metrics, picard_isize], "bam/multiqc_bam.html")
def multiqc_bam (infiles, outfile):
    statement = '''multiqc -n multiqc_bam.html -o bam bam > %(outfile)s.log 2>&1'''
    P.run(statement)

@follows (mkdir("counts"))
@merge (hisat2, "counts/gene.counts")
def feature_counts (infiles, outfile):
    infile_list = " ".join(infiles)
    statement = '''featureCounts %(feature_counts_options)s -T %(feature_counts_threads)s -a %(feature_counts_gtf)s -o %(outfile)s %(infile_list)s > %(outfile)s.log 2>&1'''
    P.run(statement, job_threads = PARAMS["feature_counts_threads"])

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
