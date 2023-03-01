README for RSVA Illumina Processing

overview 
- two step process to make and filter consensus sequence
1. single .sh file that has several components
	-- copies fastq files and removes lane info from file names (change_names_nextseq.py)
	-- creates and submits a series of bash scripts (one for each sample) that trims adapters and primers, and uses IRMA to generate consensus sequences (IRMA_RSVA_batch_script.py)
2. python script that concatenates and filters consensus sequence on completeness (filter_consensus.py)

Requirements
- IRMA path needs to be in bash profile
	-- export PATH=$PATH:/nfs/turbo/umms-alauring/shared_projects/IAV_illumina/flu-amd
- Need to have biopython and pandas installed in python/3.9.1

Steps
1. make a folder for each run inside RSVA_Illumina
2. open (vi or nano) IRMA.sh and change run name info-  ("Sequencing plate")
	- lines below are the ones that need changing
		cd Sequencing_Plate
		-s  nfs/turbo/umms-alauring/raw_data/2022/Sequencing_Plate 
		--prefix Sequencing_Plate
3. sh IRMA.sh
4. Wait for all jobs to finish (there will be 96 jobs for a full plate)
5. Change --prefix to correct run name and copy and paste into the terminal

"""
python IRMA_filter_consensus_coverage.py \

    --prefix 20210712_Nanopore_Run_33 \
    
    --mode test \
   
    --min 13763

"""
