# bioinf-notes
Creating a file with useful commands in Bioinformatics, some commands are my own and some are taken from various sources
This file is a reference only

## Contents
- [Sources](#sources)
- [RNASeq pipeline](#rnaseq-pipeline)
- [miRNAseq Pipeline](#mirnaseq-pipeline)
- [Redirect output of a command to a file](#redirect-output-of-a-command-to-a-file)
- [samtools](#samtools)
- [parsing gencode GTF file](#parsing-gencode-gtf-file)

## Sources
* <http://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file>
* <https://www.biostars.org/p/56246/>
* Pevsner, Jonathan. Bioinformatics And Functional Genomics. Hoboken: John Wiley and Sons, 2015. Print.
* <http://bib.oxfordjournals.org/content/early/2015/04/17/bib.bbv019.full>
* <https://wikis.utexas.edu/display/bioiteam/Alternative+Applications+of+RNA-seq>

## RNASeq pipeline (for Ion Torrent reads)
[[back to top](#contents)]

remove adaptor sequences

    cutadapt -m 16 -b GGCCAAGGCG -o adaptorTrim.fastq input.fastq
 
Run STAR on long reads

    STAR --genomeDir <path to STAR index> --runThreadN 12 --readFilesIn <path to fq files> --outSAMunmapped unmappedSTAR.sam --outReadsUnmapped Fastx --chimSegmentMin 18 --chimScoreMin 12

Sort SAM for STAR aligned reads

    samtools view -bS Aligned.out.sam | samtools sort - -o alignedSTAR.bam

Run bowtie2 aligner on STAR unmapped reads


    bowtie2 --local --very-sensitive-local -p 8 -q --mm -x /results/plugins/scratch/RNASeqAnalysis/hg19/bowtie2/bowtie2 -U      Unmapped.out.mate1 --un sbt2_unmap.fq  | samtools view -bS - | samtools sort - unmapped_remapBowtie2

Merge and index STAR and bowtie2 aligned reads

    java -jar picard.jar MergeSamFiles USE_THREADING=true MSD=true AS=true I=alignedSTAR.bam I=unmapped_remapBowtie2.bam O=....STARBowtie2.bam

map unmapped reads to rRNA using bowtie2

    bowtie2 --local --very-sensitive-local -p 8 -q --mm -x <path to bowtie2 index for rRNA> -U sbt2_unmap.fq | samtools view -bS - | samtools sort - -o rRNA.bam
    
*To identify human rRNA
RefSeq sequences from GenBank, follow the following steps. (1) From the home
page of NCBI, navigate to NCBI Nucleotide and restrict the search to human using
the search builder. (2) Currently (May 2015) there are nearly 11 million entries.
Click “rRNA” under the “molecule types” filter. (3) There are now 30 RefSeq
entries corresponding to 5.8S rRNA (e.g., NR_003285; 156 base pairs), 28S rRNA
(NR_003287; 5070 base pairs), 18S rRNA (NR_003286; 1869 base pairs), and 45S
rRNA (NR_046235; 13,357 base pairs). For each, the chromosomal assignment is to
the acrocentric p-arms.*


get Alignment summary metrics using picard
    
    java -jar picard.jar CollectAlignmentSummaryMetrics I=STARBowtie2.bam O=STARBowtie2.alignmentSummary.txt R=hg19.fasta LEVEL=ALL_READS

get RNASeq metrics using picard

    java -jar picard.jar CollectRnaSeqMetrics REF_FLAT=refFlat RIBOSOMAL_INTERVALS=rRNA.interval STRAND=FIRST_READ_TRANSCRIPTION_STRAND MINIMUM_LENGTH=100 LEVEL=ALL_READS I=....STARBowtie2.bam R=hg19.fasta O=...STARBowtie2.RNAmetrics.txt

gene counts using HTSeq

    samtools view -F 4 STARBowtie2.bam | htseq-count -q -t exon -i gene_name -m intersection-nonempty - /annotations/hg19/gene.gtf  >  STARBowtie2.gene.count

  - To combine multiple count files, use join
                       
        join count1.gff count2.gff| join - count3.gff | join - count4.gff |join - count5.gff|join - count6.gff > gene_counts_HTseq.gff

Run cufflink, cuffmerge and cuffdiff on aligned reads

    cufflinks -q -p 12 -m 100 -s 60 -G /annotations/hg19/gene.gtf -M /annotations/hg19/rRNA_mask.gtf   --library-type fr-secondstrand --max-bundle-length 3500000   -o output_cufflinks --no-update-check ....STARBowtie2.bam

## miRNAseq pipeline
[[back to top](#contents)]

- Optimal parameters for different aligners <http://bib.oxfordjournals.org/content/early/2015/04/17/bib.bbv019.full>

        BWA 0.7.4: bwa aln -n 1 -o 0 -e 0 -k 1 -t 4

        BWA 0.7.4 (0 mismatch in seed): bwa aln -n 1 -o 0 -e 0 -l 8 -k 0 -t 4
        
        Bowtie 0.12.9: bowtie -n 1 -l 8 -a --best --strata --phred33- quals
        
        Bowtie 0.12.9 (0 mismatch in seed): bowtie -n 0 -l 8 -a --best --strata --phred33-quals
        
        Bowtie2 2.1.0: bowtie2 --local -p 8 -q --phred33 -D 20 -R 3 -N 0 -L 8 -i S,1,0.50
        
        Novoalign 3.00.05: novoalign -a TGGAATTCTCGGGT GCCA AGG -l 15 -t 30 -r A


Build reference with bowtie

  - Get reference fasta from mirbase (then filter for organism) or should be in folder for organism downloaded from iGenome collection

        bowtie2-build mature_hsa.fa mature_hsa
    
Bowtie2 Local Alignment

  - We would like to use an alignment strategy that can intelligently ignore the parts that won't align to a reference (the 'adapter') and align correctly the parts that align well.  This is called a 'local' alignment, in contrast to a 'global' alignment, which would count the 'mismatches' in the adapter against the alignment score. Bowtie2 is a local-alignment-capable aligner 
  
            bowtie2 --local -N 1 -L 16 -x mirbase/mature_hsa -U human_mirnaseq.fastq.gz -S human_mirnaseq.sam
            
miRNA Profiling with SAMTools

  - Because each contig in our reference is a feature, we do not need a true annotation to quantify our miRNAs.  Instead, we just need to count how many reads aligned to each contig, and sort them accordingly. Samtools has a simple utility called idxstats that reports exactly this.  The following commands will produce this information by converting the SAM file to BAM, sorting and indexing it, then running idxstats: 
  
            samtools view -bS human_mirnaseq.sam > human_mirnaseq.bam
            samtools sort human_mirnaseq.bam human_mirnaseq.sorted
            samtools index human_mirnaseq.sorted.bam
            samtools idxstats human_mirnaseq.sorted.bam > idxstats.txt

  - look at the top ten most abundant miRNAs
  
            sort -r -n -k 3 idxstats.txt | head


## Redirect output of a command to a file
[[back to top](#contents)]

The standard error stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, it gets overwritten.

    command 2> output.txt
    
Both the standard output and standard error stream will be redirected to the file only, nothing will be visible in the terminal. If the file already exists, it gets overwritten.

    command &> output.txt
    
The standard output stream will be copied to the file, it will still be visible in the terminal. If the file already exists, it gets overwritten.

    command | tee output.txt
    
The standard output stream will be copied to the file only, it will still be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

    command | tee -a output.txt
    
Example with bowtie2 (saving alignment stats in log file)
    
    bowtie2 --local -p 8 -x genomePrefix -U reads.fq --un unmapped.fq -S aligned.sam 2>bowtie2.log 

## samtools
[[back to top](#contents)]

Get unmapped reads from a bam file (output in sam)

    samtools view -f 4 file.bam > unmapped.sam

Get unmapped reads from a bam file (output in bam)

    samtools view -b -f 4 file.bam > unmapped.bam

Get only mapped reads (use the parameter 'F', which works like -v of grep and skips the alignments for a specific flag)

    samtools view -b -F 4 file.bam > mapped.bam     

Get the unique reads (a single read mapping at one best position)

    samtools view -b -q 1 file.bam > unique.bam 

- Method is debatable (use MAPQ of 5 or 10...explained below)
    - *From Devon Ryan in biostars post <https://www.biostars.org/p/101533/>
'Bowtie2 will give an alignment a MAPQ score of 0 or 1 if it can map equally well to more than one location. Further, there's not always a perfect correspondence between the MAPQ of a read and the summary metrics it prints at the end (I'd need to go through the code again to determine how it arrives at the printed summary metrics, that's not documented anywhere). Finally, you would be well served to completely abandon the concept of "uniquely mapped". It is never useful and is always misleading, since you're simply lying to by labeling something unique. You're better served by simply filtering on a meaningful MAPQ (5 or 10 are often reasonable choices), which has the benefit of actually doing what you want, namely filtering according the the likelihood that an alignment is correct.'*

## parsing gencode GTF file
[[back to top](#contents)]

- <https://www.gencodegenes.org>

Get all "gene" lines

    awk '{if($3=="gene"){print $0}}' gencode.gtf
    
Get all "protein-coding transcript" lines: 

    awk '{if($3=="transcript" && $20=="\"protein_coding\";"){print $0}}' gencode.gtf
    
Get level 1 & 2 annotation (manually annotated) only:
    
    awk '{if($0~"level (1|2);"){print $0}}' gencode.gtf
    
Parse using perl <https://www.gencodegenes.org/data_format.html>

Get all miRNA gene names:

    awk '{if($20 == "\"miRNA\";"){print $0}}' gencode.v19.annotation.gtf | cut -d ";" -f 5 - | awk -F " " '{print $2}' - | sort | uniq > miRNA.genes
