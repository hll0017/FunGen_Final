#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft formate
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript level
## HiSat2  Indexing   InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort          Input: Alignment files,  .sam
##                                                  Output: Sorted  .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##    queue: class or medium
##    core: 6
##    time limit (HH:MM:SS): 04:00:00 
##    Memory: 12gb
##    
###############################################

#### Load all the programs you are going to use in this script.
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread
#module load gffcompare


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## Replace the [#] with paths to define these variable
MyID=aubclsd0322          ## Example: MyID=aubtss

WD=/scratch/$MyID/RNAseq            ## Example:/scratch/$MyID/PracticeRNAseq  
CD=/scratch/$MyID/RNAseq/CleanData            ## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired files are located from the last script
REFD=/scratch/$MyID/RNAseq/DaphniaRefGenome          ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains the indexed reference genome for the garter snake
MAPD=/scratch/$MyID/RNAseq/Map_HiSat2           ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/scratch/$MyID/RNAseq/Counts_StringTie       ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/RNAseq/Counts_H_S_2025      ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S

REF=GCF_000002285.5_Dog10K_Boxer_Tasha_genomic                  ## This is what the "easy name" will be for the genome reference

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD

### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
cp /home/$MyID/FunGen_Final/ncbi_dataset/data/GCF_000002285.5/$REF.fna .
cp /home/$MyID/FunGen_Final/ncbi_dataset/data/GCF_000002285.5/genomic.gtf .

###  Identify exons and splice sites on the reference genome
# gffread genomic.gff -T -o ${REF}.gtf               ## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
hisat2_extract_splice_sites.py genomic.gtf > ${REF}.ss
hisat2_extract_exons.py genomic.gtf > ${REF}.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Dog10K_Boxer_Tasha_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

## Move to the directory for mapping
cd ${MAPD}

## Copy SRR_IDs.txt to Mapping Directory
cp /home/$MyID/FunGen_Final/RNAseq_Samples/SRR_IDs.txt .

## Create list of fastq files to map.
SRR_IDs=$(cat SRR_IDs.txt)

## process the samples in the list, one by one using a while loop
for SRR in ${SRR_IDs[@]}
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 12 --dta --phred33       \
    -x "${REFD}"/Dog10K_Boxer_Tasha_index       \
    -U "${CD}"/"$SRR".fastq \
    -S "$SRR".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 12 -bS "$SRR".sam > "$SRR".bam

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 12  "$SRR".bam  -o  "$SRR"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat "$SRR"_sorted.bam   > "$SRR"_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 

	mkdir "${COUNTSD}"/"$SRR"
	stringtie -p 12 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$SRR"/"$SRR".gtf -l "$SRR"   "${MAPD}"/"$SRR"_sorted.bam

done

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt ${RESULTSD}

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix. 
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/FunGen_Final/prepDE.py3 .

prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv ${RESULTSD}
## move these results files to your personal computer for downstream statistical analyses in R studio.
