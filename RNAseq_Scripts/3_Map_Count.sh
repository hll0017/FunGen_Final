#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to index the reference genome using HiSat2, map SRR reads to the referenced genome with HiSat2,
##	convert .sam to .bam using SAMtools, count reads and create a count matrix using StringTie and prepDE.py
## Required files: SRR.fastq files, SRR_IDs.txt, GCF_000002285.5_Dog10K_Boxer_Tasha_genomic reference genome, prepDE.py script
## HiSat2 : https://daehwankimlab.github.io/hisat2/
## SAMtools : https://www.htslib.org/
## StringTie : https://ccb.jhu.edu/software/stringtie/
## Recommended
##    core: 12
##    time limit (HH:MM:SS): 08:00:00 
##    Memory: 240gb
###############################################

#### Load modules
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
## Make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsd0322

## Create variables for working directory (WD), clean data directory (CD), reference genome directory (REFD), mapping directory (MAPD), counts directory (COUNTSD)
## Results directory (RESULTSD), and reference genome (REF)
WD=/scratch/$MyID/RNAseq
CD=/scratch/$MyID/RNAseq/CleanData
REFD=/scratch/$MyID/RNAseq/DogReferenceGenome
MAPD=/scratch/$MyID/RNAseq/Map_HiSat2
COUNTSD=/scratch/$MyID/RNAseq/Counts_StringTie
RESULTSD=/home/$MyID/RNAseq/Counts_H_S_2025

REF=GCF_000002285.5_Dog10K_Boxer_Tasha_genomic

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD

### Copy the reference genome (.fasta) and the annotation file (.gff3) to this REFD directory
cp /home/$MyID/FunGen_Final/ncbi_dataset/data/GCF_000002285.5/$REF.fna .
cp /home/$MyID/FunGen_Final/ncbi_dataset/data/GCF_000002285.5/genomic.gff ./$REF.gff

##  Identify exons and splice sites on the reference genome
## gffread converts the annotation file from .gff3 to .gft formate for HiSat2 to use.
gffread $REF.gff -T -o $REF.gtf
hisat2_extract_splice_sites.py $REF.gtf > $REF.ss
hisat2_extract_exons.py $REF.gtf > $REF.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss $REF.ss --exon $REF.exon $REF.fna ${REF}_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

## Move to the directory for mapping
cd $MAPD

## Copy SRR_IDs.txt to Mapping Directory
cp /home/$MyID/FunGen_Final/RNAseq_Samples/SRR_IDs.txt .

## Create list of fastq files to map.
SRR_IDs=$(cat SRR_IDs.txt)

## process the samples in the list, one by one using a for loop
for SRR in ${SRR_IDs[@]}
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 12 --dta --phred33 -x $REFD/${REF}_index -U $CD/$SRR.fastq -S $SRR.sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 12 -bS $SRR.sam > $SRR.bam

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 12  $SRR.bam -o ${SRR}_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat ${SRR}_sorted.bam   > ${SRR}_Stats.txt

  ### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model.
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 

  mkdir "${COUNTSD}"/"$SRR"
  stringtie -p 12 -e -B -G  "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$SRR"/"$SRR".gtf -l "$SRR"  "${MAPD}"/"$SRR"_sorted.bam

done

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt $RESULTSD

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix. 
 ## Move to the counts directory
cd $COUNTSD
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/$MyID/FunGen_Final/RNAseq_Scripts/prepDE.py3 .

prepDE.py3 -i $COUNTSD

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv $RESULTSD
## move these results files to your personal computer for downstream statistical analyses in R studio.
