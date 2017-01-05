# bioinf-notes
Creating a file with useful commands in Bioinformatics, some commands are my own and some are taken from various sources
This file is a reference only

## Contents
- [Sources](#sources)
- [RNASeq pipeline](#rnaseq-pipeline)
- [TopHat and CuffLinks Sample Protocol](#tophat-and-cufflinks-sample-protocol)
- [miRNAseq Pipeline](#mirnaseq-pipeline)
- [Using featureCounts from subread package](#using-featurecounts-from-subread-package)
- [samtools](#samtools)
- [parsing gencode GTF file and examining GTF files](#parsing-gencode-gtf-file-and-examining-gtf-files)
- [Redirect output of a command to a file](#redirect-output-of-a-command-to-a-file)
- [Extract file name in unix loops](#extract-file-name-in-unix-loops)
- [Automate backup in Linux with cron and rsync](#automate-backup-in-linux-with-cron-and-rsync)
- [Set up NFS server](#set-up-nfs-server)
- [get files from ftp server or http using wget rsync mget](#get-files-from-ftp-server-or-http-using-wget-rsync-mget)
- [download raw sequence data from GEO-SRA](#download-raw-sequence-data-from-GEO-SRA)


## Sources
* <http://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file>
* <https://www.biostars.org/p/56246/>
* Pevsner, Jonathan. Bioinformatics And Functional Genomics. Hoboken: John Wiley and Sons, 2015. Print.
* <http://bib.oxfordjournals.org/content/early/2015/04/17/bib.bbv019.full>
* <https://wikis.utexas.edu/display/bioiteam/Alternative+Applications+of+RNA-seq>
* <http://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash>
* <https://www.marksanborn.net/howto/use-rsync-for-daily-weekly-and-full-monthly-backups/>
* <http://bioinf.wehi.edu.au/featureCounts/>
* <https://help.ubuntu.com/community/SettingUpNFSHowTo>
* <http://www.genomicscode.org/2012/09/how-to-quickly-download-illumina.html>
* <https://www.biostars.org/p/111040/>

## RNASeq pipeline
[[back to top](#contents)]

*Optimized for Ion Torrent Reads*

Obtain genome and annotations from Illumina iGenomes through ftp
    
    ftp
    open ussd-ftp.illumina.com
    *Name*: igenome
    *password*: G3nom3s4u
    cd Homo_sapiens/UCSC/hg38/
    get Homo_sapiens_UCSC_hg38.tar.gz
    
    Can also use Filezilla (faster)
    
    *Host:* ussd-ftp.illumina.com
    *Username:* igenome
    *Password:* G3nom3s4u



Create link to index files and gtf gile (example code)

    ln –s ~/data/hg19/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome* . # the * specifies all files in that directory beginning with genome. The final . symbol indicates that the link should be placed here in the current directory

    
    ln –s ~/data/hg19/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf .


  - Or export path in .bashrc or .zshrc
  
        export STARINDEX="/References_data/References_genome/Homo_sapiens/UCSC/hg19/Sequence/STAR"
        
        export GENCODE="/References_data/References_genome/Homo_sapiens/gencode.v19.annotation.gtf"
        
  - Directory or file can then be called with **$STAR** and **$GENCODE** in command line

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


## TopHat and CuffLinks Sample Protocol
[[back to top](#contents)]

1. Obtain genome and gtf files from iGenome or Ensemble

        cp ~/Downloads/Drosophila_melanogaster_Ensembl_BDGP5.25.tar ~/data/
        tar xopf Drosophila_melanogaster_Ensembl_BDGP5.25.tar
        ln –s ~/data/Drosophila_melanogaster/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome* 
        ln –s ~/data/Drosophila_melanogaster/Ensembl/BDGP5.25/Annotation/Genes/genes.gtf .
        head -2 genes.gtf
        grep –c '@' GSM794483_C1_R1_1.fq
        
2. TopHat to Map Reads to a Reference Genome

        tophat -p 8 -G genes.gtf -o C1_R1_thout genome C1_R1_1.fq C1_R1_2.fq
        samtools flagstat accepted_hits.bam

3. Cufflinks to Assemble Transcripts

        cufflinks -p 8 -o C1_R1_clout C1_R1_thout/accepted_hits.bam
        
        
    The Cufflinks outputs are sent to a folder with text files listing the loci and lengths of
    genes, transcripts, and isoforms. Next, create a file called assemblies.txt. This lists the assembly file for each sample.
    
        find . -name transcripts.gtf > assembly_list.txt
        
4. Run Cuffmerge on all the assemblies. 

    This generates a single merged transcriptome annotation. The –s option specifies the genomic DNA sequences for the reference, while –g genes.gtf is an optional reference GTF file.
    
        cuffmerge -g genes.gtf -s genome.fa -p 8 assembly_list.txt

5. Cuffdiff to Determine Differential Expression

        cuffdiff -o diff_out -b genome.fa -p 8 -L C1,C2 -u merged_asm/merged.gtf ./C1_R1_thout/accepted_hits.bam,./C1_R2_thout/accepted_hits.bam,./C1_R3_thout/accepted_hits.bam ./C2_R1_thout/accepted_hits.bam,./C2_R3_thout/accepted_hits.bam,./C2_R2_thout/accepted_hits.bam
        
6. CummeRbund to Visualize RNA-seq Results in R
        
        Rstudio
        source("http://bioconductor.org/biocLite.R")
        biocLite("cummeRbund")
        library(cummeRbund)
        
        cuff_data <- readCufflinks('diff_out')



## miRNAseq pipeline
[[back to top](#contents)]

- Basic Parameters to consider

    - Alignment parameters 

     - Minimum percent identity (%) : 96
     - Mismatches allowed (0-2) : 1

    - Output Parameters 

     - Minimum match length (bp) : 10
     - Max number of matches to be reported per read : 5

     - Ignore reads with more than _ valid matches : 5

      - Trimming parameters

        - Fixed Trimming
        
          - Number of bases to trim from 3' end : 0
          - Number of bases to trim from 5' end : 0

        - Quality Trimming 

          - Trim 3'end with average quality less than : 10
          
       - Perform Adapter Trimming : No
       - Trim poorly aligned portion at 3'end : No

    - Screen against Ribosomal RNA contaminants


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

Consider Subread package for miRNA seq alignment too as they have parameters optimized for smallRNA seq. Check users manual

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

Get only paired-end reads from bam file

    samtools view -bf 1 foo.bam > foo.paired-end.bam
    
Get only single-end reads from bam file

    samtools view -bF 1 foo.bam > foo.single-end.bam

Examine a few lines of BAM alignment file.

    samtools view -x accepted_hits.bam | less

- Spliced sequences

  - The 6th BAM file field is the CIGAR string which tells you how your query sequence mapped to the reference.
   
  - The CIGAR string "58M76N17M" represents a spliced sequence. The codes mean:
    - 56M-the first 58 bases match the reference
    - 76N-there are then 76 bases on the reference with no corresponding bases in the sequence (an intron)
    - 17M-the last 17 bases match the reference
   
  - Count spliced sequences
   
            samtools view accepted_hits.bam | cut -f 6 | grep 'N' | wc -l


## Using featureCounts from subread package
[[back to top](#contents)]

**featureCounts is a way more faster than HTSeq-counts**

Summarize a *single-end* read dataset using 5 threads:

        featureCounts -T 5 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.sam
        
Summarize a BAM format dataset:

        featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam

Summarize multiple datasets at the same time:

        featureCounts -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam

Perform strand-specific read counting (use '-s 2' if reversely stranded):
        
        featureCounts -s 1 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_SE.bam

Summarize *paired-end* reads and count fragments (instead of reads):

  - From subread manual: Reads may be paired or unpaired. If paired reads are used, then each pair of reads defines
a DNA or RNA fragment bookended by the two reads. In this case, featureCounts can be
instructed to count fragments rather than reads. featureCounts automatically sorts reads by
name if paired reads are not in consecutive positions in the SAM or BAM file, with minimal
cost. Users do not need to sort their paired reads before providing them to featureCounts.

        
        featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

Summarize multiple paired-end datasets:
        
        featureCounts -p -t exon -g gene_id -a annotation.gtf -o counts.txt library1.bam library2.bam library3.bam

Count the fragments that have fragment length between 50bp and 600bp only:
        
        featureCounts -p -P -d 50 -D 600 -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

Count those fragments that have both ends mapped only:
        
        featureCounts -p -B -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

Exclude chimeric fragments from fragment counting:
        
        featureCounts -p -C -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

## parsing gencode GTF file and examining GTF files
[[back to top](#contents)]

- <https://www.gencodegenes.org>

Get all "gene" lines

    awk '{if($3=="gene"){print $0}}' gencode.gtf
    
Get all "protein-coding transcript" lines: 

    awk '{if($3=="transcript" && $20=="\"protein_coding\";"){print $0}}' gencode.gtf
    
Get level 1 & 2 annotation (manually annotated) only:
    
    awk '{if($0~"level (1|2);"){print $0}}' gencode.gtf
    
Parse using perl <https://www.gencodegenes.org/data_format.html>

Get gene annotation:

    awk '{if($3=="gene"){print $0}}' gencode.gtf | cut -f 9 | cut -d ' ' -f 2,6,10 >  gencode_annotation.txt 
    
Get all miRNA gene names:

    awk '{if($20 == "\"miRNA\";"){print $0}}' gencode.v19.annotation.gtf | cut -d ";" -f 5 - | awk -F " " '{print $2}' - | sort | uniq > miRNA.genes
    
If using genes.gtf from Gencode.genes folder from hg38 iGenome and want a text file of gene_id and gene_name

    awk '{if($9=="gene_id"){print $0}}' genes.gtf | cut -f 9 | cut -d ";" -f 1,2 | cut -d " " -f 2,4 > Genes.gencode_gene2symbol.txt
    
Other
    
    less $BI/ngs_course/tophat_cufflinks/reference/genes.gtf  # :q to exit
    cat $BI/ngs_course/tophat_cufflinks/reference/genes.gtf | head
    cat $BI/ngs_course/tophat_cufflinks/reference/genes.gtf | cut -f 1-8 | more
    cat $BI/ngs_course/tophat_cufflinks/reference/genes.gtf | cut -f 9 | more

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


## Extract file name in unix loops
[[back to top](#contents)]

*FILE="example.tar.gz"*

Extract 'example'

    echo "${FILE%%.*}"
    
Extract 'example.tar'

    echo "${FILE%.*}"
    
Extract 'tar.gz'

    echo "${FILE#*.}"
        
Extract 'gz'        
        
    echo "${FILE##*.}"
    
Trim the shortest match from the end
    
    ${variable%pattern}

Trim the longest match from the beginning

    ${variable##pattern}

Trim the longest match from the end

    ${variable%%pattern}

Trim the shortest match from the beginning
    
    ${variable#pattern}

## Automate backup in Linux with cron and rsync
[[back to top](#contents)]

Open chrontab text file with this command

        crontab -e
        
Add lines of code at the bottom (examples below)

    30 17 * * * rsync –av /path/to/source /home/mark/rsync/daily
    00 18 * * 5 rsync –av --delete /home/mark/rsync/daily /home/mark/rsync/weekly
    00 6 1 * * tar -cvjf /home/mark/rsync/monthly/monthly_$(date +%Y%m%d).tar.bz2 /home/rsync/daily/
    
In the above example cron setup will backup daily at 5:30PM, Backup every Friday at 6:00PM and Do the full backup on the first of each month at 6:00AM

Adding the *--delete* flag: By default rsync does not delete files in the backup if they were deleted in the source. This is rsync’s way of protecting you from deleting the entire backup directly on accident.

Cron by default sends emails with the output of the command. If you don’t want to get emails you can pipe the cron comands to /dev/null or to a log file

    59 */6 * * * 30 17 * * * rsync –av --delete /path/to/source /home/mark/rsync/daily >> log.file
    

## Set up NFS server
[[back to top](#contents)]

### On server end

    sudo apt-get install nfs-kernel-server 
    gedit /etc/exports
    
To export our directories to a local network 192.168.1.0/24, we add the following line to /etc/exports

    /home            192.168.0.1/24(rw,sync,no_root_squash,no_subtree_check)
    
*Note: Note that when locking down which clients can map an export by setting the IP address, you can either specify an address range (as shown above) using a subnet mask, or you can list a single IP address followed by the options. Using a subnet mask for single client's full IP address is **not** required. Just use something like 192.168.1.123(rw). There are a couple options for specifying the subnet mask. One style is 255.255.255.0. The other style is /24 as shown. Both styles should work. The subnet mask marks which part of IP address must be evaluated. *

Export file 

    exportfs -a

Restart the service 

    sudo service nfs-kernel-server restart
    sudo service idmapd restart  #Optional
    sudo modprobe nfs  #if error message 'mount.nfs4: No such device'
    
### On client end

    sudo apt-get install nfs-common rpcbind
    
create the NFS directory mount point    

    sudo mkdir -p /mnt/nfs/home
    
mount the NFS shared content in the client machine

    sudo mount -t nfs4 192.168.0.1/home /mnt/nfs/home/
    
Cross check 
    
    sudo mount -t nfs
    showmount -e <NFS server name>
    
Restart nfs-common (optional)

    sudo service idmapd restart 
    
Unmount

    sudo umount /mnt/nfs/home
    
For unmount if umount is hanging

    sudo umount -f -l /mnt/nfs/home
    sudo service nfs-kernel-server restart (on *server* side, optional)
    sudo service idmapd restart (on *client* side, optional)
 
- Helpful links

  - <https://help.ubuntu.com/14.04/serverguide/network-file-system.html>  
  - <https://help.ubuntu.com/community/SettingUpNFSHowTo#NFS_quick_start>  
  - <https://help.ubuntu.com/community/NFSv4Howto>  
  - <https://www.howtoforge.com/nfs-server-on-ubuntu-14.10>  
  - <https://www.digitalocean.com/community/tutorials/how-to-set-up-an-nfs-mount-on-ubuntu-14-04>  
  - <http://unix.stackexchange.com/questions/106122/mount-nfs-access-denied-by-server-while-mounting-on-ubuntu-machines>  
  - <http://stackoverflow.com/questions/40317/force-unmount-of-nfs-mounted-directory>


## get files from ftp server or http using wget rsync mget

    wget -m ftp://caftpd.nci.nih.gov/pub/dcc_target/RT/miRNA-seq/L3/expression/BCCA/
    
Using an rsync command to download the entire directory:

    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ .
    
Using rsync for a single file, e.g. gc5Base.txt.gz

    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gc5Base.txt.gz .
    
wget, all files
    
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/*'
        
wget, single file

    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gc5Base.txt.gz' 
        -O gc5Base.txt.gz
        
mget        

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory) 
    
## download raw sequence data from GEO SRA

https://www.ncbi.nlm.nih.gov/sra/
https://www.ncbi.nlm.nih.gov/geo/
http://www.ebi.ac.uk/ena/

After finding samples on SRA or GEO, find *SRR* number for sample on the NCBI SRA website, then use SRA toolkit on cmd line
e.g. for "SRX1210960: GSM1872102: Panc10.05 #1; Homo sapiens; RNA-Seq", link https://www.ncbi.nlm.nih.gov/sra/?term=SRX1210960

    prefetch -v SRR2313114
    
Note where the sra file is downloaded (by default to /home/[USER]/ncbi/public/sra/.) and then convert to fastq with something like the following.

    fastq-dump --outdir /opt/fastq/ --split-files /home/[USER]/ncbi/public/sra/SRR2313114.sra

    
If you just want to download X number of raw (fastq) reads to standard output from a particular run you can use a command like the following. This can be useful to just take a quick look at some reads, or obtain some reads for testing purposes or just check whether the SRA toolkit is even working for you

    fastq-dump -X 5 -Z SRR925811
