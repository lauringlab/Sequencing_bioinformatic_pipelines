#!/bin/bash

#SBATCH --job-name=IBV_IRMA
#SBATCH --mail-type=BEGIN,END
#SBATCH --mem-per-cpu=15000m 
#SBATCH --time=1-10:10:00
#SBATCH --account=alauring99
#SBATCH --partition=standard


# change any where it says run_name to the actual plate name (e.g 20220510_IAV_Illumina_Run_20)
# change name_of_sequencing_run to folder name that contains the raw sequencing data (e.g 20220523_IAV_Illumina_Run20_21)
# change barcode map to name of csv file (20220510_IAV_Illumina_Run_20_BarcodeMap.csv)
	## should have 2 columns (Processing_Plate and sample_name) 
	## All sample names with _ need to be changed to -
		# NP20_IAV_HeLa_S3_RNACtrl to NP20-IAV-HeLa-S3-RNACtrl
# need to have IRMA path in bash_profile
## EXAMPLE
	#module load python/3.9.1
	#cd run_name
	#python ../cp_fastq_files.py --sequencing_run name_of_sequencing_run --map barcode_map --mode run
	#python ../change_names_miseq.py -s data/fastq -f data/fastq_renamed -run
	#python ../Run_IRMA_Pipeline_IAV_illumina.py --prefix run_name --mode run --min-length 11819
# to submit save file and then command sbatch IRMA.sh


module load python/3.9.12
module load Bioinformatics 
module load cutadapt
module load R/4.2.0
cd test
mkdir data/cutadapt
#python ../cp_fastq_files.py --sequencing_run 20220926_IAV_Illumina_Run_28 --map 20220926_IAV_Illumina_Run_28.csv --mode run
#python ../change_names_miseq.py -s data/fastq -f data/fastq_renamed -run
python ../Run_IRMA_Pipeline_IBV_illumina.py --prefix test --mode run --min-length 11819
