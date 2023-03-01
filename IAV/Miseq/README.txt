READ ME for IAV Miseq illumina processing

Overview
- A single shell script ( IRMA_miseq.sh) calls several python scripts to do the following:
	- copies relevant fastq files from sequencing run (2+ sequencing plates) to sequencing plate folder ("cp_fastq_files.py")
		- I refer to sequencing run as all pooled samples that share a single miseq run
		- sequencing plate refers to the 96 well plate that RNA harvest and library prep occurs on
	- renames fastq files to get rid of sequencing info in file name ("change_names_miseq.py")
	- uses cutadapt to trim sequencing adapters off ends of reads, uses IRMA to make consensus sequence, organizes files for further analysis ("Run_IRMA_Pipeline_IAV_illumina.py")
		- uses flu-amd/IRMA_RES/modules/FLU/config/FLU.sh as a configuration file 
			- make any changes to IRMA's configuration in this file (ex. required PHRED score)
			- all default configurations are in flu-amd/IRMA_RES/modules/FLU/init.sh
	-
Requirements
 - need to make a folder for each sequencing plate
 - need a barcode map (Run_name_BarcodeMap.csv) within the above folder
	- two columns (Processing_Plate,sample_name)
	- replace any "_" in sample names with "-"
		- ex: NP21_IAV_HeLa_S3_RNACtrl  to NP21-IAV-HeLa-S3-RNACtr
	- see example_BarcodeMap.csv for example barcode map
- update IRMA.sh with relevant names
	- add name of barcode map to "--map 'barcode_map' " 
	- add name of sequencing run to  "--sequencing_run 'name_of_sequencing_run' "
	- add name of sequencing plate to "cd 'run_name' " and "--prefix 'run_name' "
- before first use add IRMA to bash profile 
	- export PATH=$PATH:/nfs/turbo/umms-alauring/shared_projects/IAV_illumina/flu-amd
	
Usage: "sbatch IRMA_miseq.sh" 

