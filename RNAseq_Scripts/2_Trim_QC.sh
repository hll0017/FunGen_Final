#! /bin/bash

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data
## Required files : SRR.fastq files from 1_Download_QC.sh, SRR_IDs.txt, TruSeq3-SE.fa
## FastQC : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
## Recommended parameters
        ## core: 12
        ## time limit (HH:MM:SS): 4:00:00
        ## Memory: 64gb
###############################################


# Load modules
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

## Make variables so the directories are automatically made in YOUR directory
MyID=aubclsd0322

# Variables: raw data directory (DD), working directory(WD), and Quality after cleaning (PCQ)
WD=/scratch/$MyID/RNAseq
DD=/scratch/$MyID/RNAseq/RawData
CD=/scratch/$MyID/RNAseq/CleanData 
PCQ=PostCleanQuality

## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir $CD
mkdir $WD/$PCQ

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd $DD

### Copy SRR_IDs.txt to obtain list of sequences to trim
cp /home/$MyID/FunGen_Final/RNAseq_Samples/SRR_IDs.txt .

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
cp /home/$MyID/FunGen_Final/TruSeq3-SE.fa .

### Create a variable for SRR IDs
SRR_IDs=$(cat SRR_IDs.txt)

### Run a for loop to Trim reads with Trimmomatic and assess quality with FastQC
for SRR in ${SRR_IDs[@]}
do

        ### Run Trimmomatic in single end (SE) mode with 12 threads using phred 33 quality score format. 
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar  \
	SE -threads 12 -phred33 "$SRR".fastq $CD/"$SRR".fastq  \
        ILLUMINACLIP:TruSeq3-SE.fa:2:35:10 HEADCROP:5 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results
	fastqc $CD/"$SRR".fastq --outdir=$WD/$PCQ

done
# This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd $WD/$PCQ

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf $PCQ.tar.gz $WD/$PCQ/*

## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate the quality of the data.
