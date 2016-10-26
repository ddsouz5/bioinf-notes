# bioinf-notes
Creating a file with useful commands in Bioinformatics, some commands are my own and some are taken from various sources
This file is a reference only

## Contents
- [Sources](#sources)
- [RNA-Seq pipeline](#RNA-Seq-pipeline)
- [Redirect output of a command to a file](#Redirect output of a command to a file)
- [samtools](#samtools)

## Sources
* <http://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file>
* <https://www.biostars.org/p/56246/>

## RNA-Seq pipeline
[[back to top](#contents)]
 
Run STAR on long reads

    STAR --genomeDir <path to STAR index> --runThreadN 12 --readFilesIn <path to fq files> --outSAMunmapped unmappedSTAR.sam --outReadsUnmapped Fastx --chimSegmentMin 18 --chimScoreMin 12

Sort SAM for STAR aligned reads

    samtools view -bS Aligned.out.sam | samtools sort - -o alignedSTAR.bam

Run bowtie2 aligner on STAR unmapped reads


    bowtie2 --local --very-sensitive-local -p 8 -q --mm -x /results/plugins/scratch/RNASeqAnalysis/hg19/bowtie2/bowtie2 -U      Unmapped.out.mate1 --un sbt2_unmap.fq  | samtools view -bS - | samtools sort - unmapped_remapBowtie2

Merging and indexing STAR and bowtie2 aligned reads

    java -jar picard.jar MergeSamFiles USE_THREADING=true MSD=true AS=true I=alignedSTAR.bam I=unmapped_remapBowtie2.bam O=....STARBowtie2.bam

map unmapped reads to rRNA

    bowtie2 --local --very-sensitive-local -p 8 -q --mm -x <path to bowtie2 index for rRNA> -U sbt2_unmap.fq | samtools view -bS - | samtools sort - -o rRNA.bam

get Alignment summary metrics
    
    java -jar picard.jar CollectAlignmentSummaryMetrics I=...STARBowtie2.bam O=....STARBowtie2.alignmentSummary.txt R=hg19.fasta LEVEL=ALL_READS

get RNASeq metrics

    java -jar picard.jar CollectRnaSeqMetrics REF_FLAT=refFlat RIBOSOMAL_INTERVALS=rRNA.interval STRAND=FIRST_READ_TRANSCRIPTION_STRAND MINIMUM_LENGTH=100 LEVEL=ALL_READS I=....STARBowtie2.bam R=hg19.fasta O=...STARBowtie2.RNAmetrics.txt

gene counts

    samtools view -F4 ....STARBowtie2.bam | htseq-count -q -t exon -i gene_name - /annotations/hg19/gene.gtf  > ...STARBowtie2.gene.count

assemble isoforms? <- check

    cufflinks -q -p 12 -m 100 -s 60 -G /annotations/hg19/gene.gtf -M /annotations/hg19/rRNA_mask.gtf   --library-type fr-secondstrand --max-bundle-length 3500000   -o output_cufflinks --no-update-check ....STARBowtie2.bam

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
    
    bowtie2 --local -p 8 -x genomePrefix -U file.fw --un unmapped.fq 2>bowtie2.log

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
*From Devon Ryan in biostars post <https://www.biostars.org/p/101533/>
'Bowtie2 will give an alignment a MAPQ score of 0 or 1 if it can map equally well to more than one location. Further, there's not always a perfect correspondence between the MAPQ of a read and the summary metrics it prints at the end (I'd need to go through the code again to determine how it arrives at the printed summary metrics, that's not documented anywhere). Finally, you would be well served to completely abandon the concept of "uniquely mapped". It is never useful and is always misleading, since you're simply lying to by labeling something unique. You're better served by simply filtering on a meaningful MAPQ (5 or 10 are often reasonable choices), which has the benefit of actually doing what you want, namely filtering according the the likelihood that an alignment is correct.'*
