## For Processing RSVB sequenced by nanopore

# Overview
	# change_names_ONT.py collects and renames the fastq files for each sample. 
	#RSVB_ONT_Snakemake filters by amplicon size, primers are trimmed by cutadapt, sequences are aligned to the GISAID reference () using minimap2, and ivar generates the consensus sequence.
	# The majority of the processing will be done on great lakes. 

# Requirements
	# need to have  biopython installed on python/3.9.12
	# need to have a conda environment for snakemake with chopper installed

# Steps 	
	# make new folder for the run. Type cd into the terminal and then drag the file from the finder window to the terminal. Then hit return.


	# Move barcode map to greatlakes. Change barcode map name and user name as needed. The barcode map should be a csv file with 2 columns  "barcode,sample_id". Barcodes that were not used should be removed from this file. 
	

scp 20231030_RSVB_Nanopore_Run_3_BarcodeMap.csv bendalle@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-alauring/shared_projects/RSVB_ONT


	#log into great lakes

ssh bendalle@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-alauring/shared_projects/RSVB_ONT

	# make run folder and add barcode map on greatlakes 

cd /nfs/turbo/umms-alauring/shared_projects/RSVB_ONT
mkdir 20231030_RSVB_Nanopore_Run_3 ## add run name
mv 20231030_RSVB_Nanopore_Run_3_BarcodeMap.csv  20231030_RSVB_Nanopore_Run_3  ## add barcode map and name of the run
cd 20231030_RSVB_Nanopore_Run_3 # move the run folder
conda activate snakemake 

	#submit job with the appropriate barcode map name and run name in this order. No flags are needed
	## example sbatch ../RSVB_ONT.sh "barcode map" "run name" 

sbatch ../RSVB_ONT.sh 20231030_RSVB_Nanopore_Run_3_BarcodeMap.csv 20231030_RSVB_Nanopore_Run_3
logout # logout of great lakes. 

	#Once you have received an email that the job is complete, bring down the processed files from greatlakes. Make sure you are in the same directory as the first step. 

scp bendalle@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-alauring/shared_projects/RSVB_ONT/20231030_RSVB_Nanopore_Run_3/20231030_RSVB_Nanopore_Run_3* .
 

