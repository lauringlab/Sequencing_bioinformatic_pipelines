README for IBV Illumina Nextseq Processing 

overview 
- two step process to make and filter consensus sequence
1. single .sh file that has several components
	-- copies relevant fastq files from sequencing run (2+ sequencing plates) to sequencing plate folder ("cp_fastq_files.py")
	-- copies fastq files and removes lane info from file names (change_names_nextseq.py)
	-- creates and submits a series of bash scripts (one for each sample) that trims adapters and primers, and uses IRMA to generate consensus sequences (IRMA_IAV_batch_script.py)

			- uses flu-amd/IRMA_RES/modules/FLU/config/FLU.sh as a configuration file 
				- make any changes to IRMA's configuration in this file (ex. required PHRED score)
				- all default configurations are in flu-amd/IRMA_RES/modules/FLU/init.sh
				
2. python script that concatenates and filters consensus sequence on completeness (IBV_filter.py)

Requirements
- IRMA path needs to be in bash profile
	-- export PATH=$PATH:/nfs/turbo/umms-alauring/shared_projects/IAV_illumina/flu-amd
- Need to have biopython and pandas installed in python/3.9.12

Steps
1. make a folder for each run inside IAV_illumina
2. need a barcode map (Run_name_BarcodeMap.csv) within the above folder
	- two columns (Processing_Plate,sample_name)
	- replace any "_" in sample names with "-"
		- ex: NP21_IAV_HeLa_S3_RNACtrl  to NP21-IAV-HeLa-S3-RNACtr
		- see example_BarcodeMap.csv for example barcode map
3. open (vi or nano) IRMA.sh and change run name info-  
	- 
	- lines below are the ones that need changing
		cd   "Plate name"
		python ../cp_fastq_files.py --sequencing_run "Sequencing run" --map "Plate map" --mode run
		python ../IRMA_IBV_batch_script.py --prefix "Plate name" --mode run

4. sh IRMA_nextseq.sh
5. Wait for all jobs to finish (there will be 96 jobs for a full plate)
6. Change --prefix to correct run name and copy and paste into the terminal

"""
python IBV_filter.py \

    --prefix 20210712_Nanopore_Run_33 \
    
    --mode test \
   
    --min 11819

"""
