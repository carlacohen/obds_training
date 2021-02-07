''' This pipeline does from peakcalling to downstream analysis

goes from bam (produced by pipeline_peakcalling_mapping.py) to bed and bigwig files

Author: Charlie-George

- input = indexed bam files e.g. sample_KO_R1.bam

- outputs
    - bam with duplicates removed
    - bam filtered by samtools
    - bam with reads that overlap blacklist regions removed
    - macs2 output files e.g narrowPeak
    - bigwig files for visulisation .bw
    - plots

'''

# import packages
import sys
from ruffus import *
from cgatcore import pipeline as P

# read in paramerters from yml file
params = P.get_parameters("peakcalling_pipeline_solution_cg.yml")


@mkdir('filtered_bams')
@transform('*.bam',
            regex(r'(.*).bam'),
            r'filtered_bams/\1.rmdup.bam')
def picard_rmdup(infile, outfile):
    ''' take in bam file and remove duplicates using picard MarkDuplicates'''

    final_memory = str(int(params['picard_memory'])+ 2)+'g'
    statement = """picard -Xmx%(picard_memory)sg
                            MarkDuplicates
                            REMOVE_DUPLICATES=true
                            I= %(infile)s
                            O=%(outfile)s
                            M=%(outfile)s.metrics
                            VALIDATION_STRINGENCY=SILENT
                            && samtools index %(outfile)s
                            """

    P.run(statement,
          job_queue=params['queue'], job_memory=final_memory)

#This function is to filter reads, umapped reads
@transform(picard_rmdup,
            regex(r'filtered_bams/(.*).rmdup.bam'),
            r'filtered_bams/\1.filter_read.bam')
def filter_read(infile, outfile):
    ''' take in bam that has had duplicates removed
        - uses samtools to filter out:
            - secondary alignments -F 0x100
            - not proper pairs   -f 0x2
            - unmapped reads -F4
            - low quality score  i.e. bowtie2 multimapping reads -q40
            - all chromosomes except chrM
    '''
    # should put this in yml file
    all_chrs = ['chr1','chr2','chr3','chr4','chr5',
                'chr6','chr7','chr8','chr9', 'chr10',
                'chr11','chr12','chr13','chr14','chr15',
                'chr16','chr17','chr18','chr19', 'chrX','chrY']

    # ' '.join(all_chrs) = makes list into a string
    chrs_to_keep = ' '.join(all_chrs)

    #-h Include the header in the output
    #-F gets rid of unmapped pairs
    #-F 0x100 gets rid of not primary alignment
    #-b is bam file - we want to output bam
    # -q is needed to filter based on MAPQ scores - >40 will remove most multimapping reads for bowite2
    # https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
    statement = """samtools view -h -b -f0x2 -F 4 -F0x100 -q 20
                    %(infile)s %(chrs_to_keep)s > %(outfile)s
                    && samtools index %(outfile)s"""

    P.run(statement,job_queue=params['queue'])


    # chr to keep are listed after the infile - happy to just list these in yml file in the exercise but if there is extra time
    # think about a way using a file to generate the list the contigs (e.g. idxstats) they could automate this - e.g write a small bit of python code if they want/have time


#Now to filter blacklisted peaks/regions
#This function will remove the regions from the params file using bedtools, also index to
@transform(filter_read,
            regex(r'filtered_bams/(.*).filter_read.bam'),
            r'filtered_bams/\1.whitelist.bam')
def blacklist_filter(infile, outfile):
    ''' remove reads in blacklist befile regions from bamfile '''

    # older versions of samtools need -abam instead -a
    statement = """ bedtools intersect -v
                                -a %(infile)s -b %(blacklist)s
                                > %(outfile)s
                                && samtools index %(outfile)s """

    P.run(statement,
          job_queue = params['queue'])



# Now for peakcalling using macs2
@follows(mkdir('peakcalling'))
@transform(blacklist_filter,
            regex(r'filtered_bams/(.*).whitelist.bam'),
            r'peakcalling/\1.macs2_peaks.narrowPeak')
def peakcalling(infile, outfile):

    outfile_prefix = P.snip(outfile,'_peaks.narrowPeak')


    # this should be parameterised more and put in yml
    # bdg makes bedgraph
    # --call-summits - calls the summit possitions
    # keep duplicates because we already removed them earlier
    # if we were to merge bam files before running peakcalling we might have got 'duplicates' that originated from different files

    # remove '.whitelist.bam' from the end of infile name
    inf = P.snip(infile,'.whitelist.bam')

    statement = """macs2 callpeak
                        --format=BAMPE
                        --treatment %(infile)s
                        --verbose=10
                        --name=%(outfile_prefix)s
                        --bdg
                        --SPMR
                        --call-summits
                        --keep-dup all
                        --gsize mm     >& %(outfile)s.log"""

    # need to specify different conda environment
    P.run(statement,
                job_queue = params['queue'],
                job_condaenv='peaktools_env')


################# PART 1 OF EXERCISE
################# PART 2 OF EXERCISE

### Make bigwig file - what are we doing
### Make one for each file
### normalise using CPM, bin size
@follows(peakcalling,mkdir('bigwigs'))
@transform(blacklist_filter,
            regex(r'filtered_bams/(.*).whitelist.bam'),
            r"bigwigs/\1.bw")
def bam2bw(infile, outfile):
    '''Make bam into bigwig file

    input: Bam - needs to be indexed
    output: bigwig in CPM

    uses deeptools to generate bigwig (.bw) file
    '''

    # this should be parameterised and put in yml
    # could use --exactScaling - this is slower but more accurate
    statement = '''bamCoverage -b %(infile)s
                        --outFileName %(outfile)s -p 12
                        --outFileFormat bigwig
                        --binSize 1
                        --normalizeUsing CPM
                        --blackListFileName  %(blacklist)s
                        --extendReads
                        --effectiveGenomeSize 2407883318 '''  #from deeptools documentation readlength 75bp

    # again need specific conda env
    P.run(statement, job_threads = 12, job_condaenv='peaktools_env')


# for motif calling want to make a master peaklist for each condition
# bedtools intersect
    # keep merge peaks from all 3 replicates for WT and  then for KO
# merge files based on sample name based on regex
# collate - there is an example to point them to on the ruffus collate page  http://www.ruffus.org.uk/decorators/collate.html
@collate(peakcalling,
            regex(r'peakcalling/(.*)_(.*)_(.*).macs2_peaks.narrowPeak'),
            r"\1_\2.merged.bed")
def merge_beds(infiles,outfile):
    '''
    merged bed files together based
    on the first and second regex group

    e.g. replicates will be merged
        e.g. 'mm_Chd7KOatac_R1', 'mm_Chd7KOatac_R3','mm_Chd7KOatac_R3' -> merged
             'mm_Chd7WTatac_R1', 'mm_Chd7WTatac_R3','mm_Chd7WTatac_R3' -> merged

    infiles  = replicates of KO or WT
    -------
    e.g. infiles = ('mm_Chd7KOatac_R1', 'mm_Chd7KOatac_R3','mm_Chd7KOatac_R3')
    outfile = a single bed file

    task - merge peaks from all three files using bedtools '''

    # infiles will be ('peakcalling/mm_Chd7KOatac_R1.trimmed.chr5.macs2_peaks.narrowPeak', 'peakcalling/mm_Chd7KOatac_R2.trimmed.chr5.macs2_peaks.narrowPeak', 'peakcalling/mm_Chd7KOatac_R3.trimmed.chr5.macs2_peaks.narrowPeak')
    # outfile will be mm_Chd7KOatac.merged.beds
    print(infiles)
    print(outfile)

    # tell them they need to cat all the files together - then direct to the bedtools merge page,
    # the instructions for sorting the file are on the page
    # make tuple of infiles into a string
    files = ' '.join(infiles)
    statement = ''' cat %(files)s | sort -k1,1 -k2,2n > %(outfile)s.sorted.bed &&
                    bedtools merge -i %(outfile)s.sorted.bed > %(outfile)s '''

    P.run(statement)

### Note we don't do motif calling here
    ## could do this later -> test that homer is available
    ## meme requires equal length fragments - do this based on summit possitions
    ## these steps are a bit too complex to take us through here and take alot of time
    ## Think about including as a seperate tutorial ..../extra for code clinics


###  use deeptools to make some plots
###  make bedfiles to give to plotting command using
###  less /ifs/mirror/annotations/mm10_ensembl81/geneset.dir/coding_gene_region.bed.gz | grep chr5 > mm10_coding_genes_region.chr5.bed
###  less /ifs/mirror/annotations/mm10_ensembl81/geneset.dir/coding_gene_region.bed.gz  > mm10_coding_genes_region.bed
###  less /ifs/mirror/annotations/mm10_ensembl81/geneset.dir/coding_gene_tss.bed.gz | grep chr5 > mm10_coding_genes_tss.chr5.bed
###  less /ifs/mirror/annotations/mm10_ensembl81/geneset.dir/coding_gene_tss.bed.gz >  mm10_coding_genes_tss.bed
###  these are all in the data directory
@follows(mkdir('plots'))
@merge(bam2bw, r"plots/matrix.gz")
def makematrix(infile, outfile):
    ''' compute coverage around each gene TSS for later plotting '''

    infiles = ' '.join(infile)
    # Make heatmap
    statement = '''
                computeMatrix reference-point --referencePoint TSS \
                -b 1000 -a 1000 \
                -R mm10_coding_genes_region.chr5.bed \
                -S %(infiles)s \
                --skipZeros \
                -o %(outfile)s \
                -p 6 \
                --outFileSortedRegions %(outfile)s.sortedregions.bed'''

    P.run(statement, job_queue=params["queue"],job_condaenv='peaktools_env',job_threads=6)


## this doesn't really show anything - except two of the samples are different
## hard to see any effect of the KO - this is because its a global down regulation of accesibility
@transform(makematrix, suffix(".gz"), "profile_plot.png")
def make_profile_plot(infile, outfile):
    ''' make line plot of coverage around TSS'''
    statement = '''plotProfile -m %(infile)s \
                    -out %(outfile)s \
                    --perGroup \
                    --plotTitle ""  \
                    --refPointLabel "TSS" \
                    -T "atac read density" \
                    -z ""
                '''
    P.run(statement, job_queue=params["queue"],job_condaenv='peaktools_env',job_threads=6)


## again no real difference in accessibility around tss of gene
## might be better to calculate and plot coverage around chd7 binding sites
@transform(makematrix, suffix(".gz"), "heatmap_plot.png")
def make_heatmap_plot(infile, outfile):
    '''make heatmap of coverage around TSS'''

    statement = '''plotHeatmap -m %(infile)s \
    -out %(outfile)s \
    --colorMap RdBu \
    --whatToShow 'heatmap and colorbar' \
    '''
    #--zMin -4 --zMax 4
    P.run(statement, job_queue=params["queue"],job_condaenv='peaktools_env',job_threads=6)



if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
