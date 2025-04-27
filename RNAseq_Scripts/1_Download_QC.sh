#! /bin/bash

######### Download Data and Quality Control ############
## Purpose: The purpose of this script is to download data from the SRA and evalulate raw data quality
## Required files : SRR_IDs.txt (located in RNAseq_Samples directory)
## Recommended parameters (on ASC):
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 36gb
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsd0322  

## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/$MyID/RNAseq/RawData   			
WD=/scratch/$MyID/RNAseq				
RDQ=RawDataQuality

##  make the directories in SCRATCH for holding the raw data
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p $DD

## move to the Data Directory
cd $DD

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
	## -F 	Defline contains only original sequence name.

# Print to stdout
vdb-config --interactive

# Copy SRR_IDs.txt to working directory
cp /home/$MyID/RNAseq_Samples/SRR_IDs.txt .

# Make a variable of SRR IDs
SRR_IDs=$(cat SRR_IDs.txt)

# Import data from NCBI
for SRR in ${SRR_IDs[@]}
do
	fastq-dump -F $SRR
done

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
mkdir $WD/$RDQ
fastqc *.fastq --outdir=$WD/$RDQ

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd $WD/$RDQ
tar cvzf $RDQ.tar.gz  $WD/$RDQ/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
