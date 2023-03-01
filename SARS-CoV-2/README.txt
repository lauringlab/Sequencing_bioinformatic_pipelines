README for SARS-CoV-2 Illumina sequencing processing


### For Single plate per lane Miseq Sequencing ###
Overview
- A series of scripts that gets a consensus file, coverage file, and variant file for SARS-CoV-2 sequenced on Illumina
- Copies sequence data from miseq output and updates fastq file names (change_names_miseq.py)
- Creates consensus sequences and coverage files for each sample (Snakemake-BWA)
	- uses BWA to align the sequences to the reference genome
	- uses ivar to trim off primers 
	- uses ivar to create a consensus sequence
	- uses samtools to create coverage file
	- uses ivar_CombineAndExport.py to concatenate and rename coverage and consensus files

Dependencies
- Reference fastq file that has been indexed by BWA (nCov_WH1_ref.fasta)
- bed file of primer locations (SARS-CoV-2_V4.1_primers_alt_overlap_sites.bed)

Requirements
- mkdir for each sequencing plate
- before using first time install conda environment for snakemake
	- install conda if you don't have it already, see https://github.com/um-dang/conda_on_the_cluster/blob/master/README.md
	- conda create --name snakemake --file /nfs/turbo/umms-alauring/shared_projects/ivyic/conda_requirements.txt
	- conda activate snakemake #(check to see if it works)
- change "Sequencing_Plate" to name of run in SARS-CoV-2_Illumina_miseq.sh

Usage: sbatch SARS-CoV-2_Illumina_miseq.sh 


###  For Multiple plates per lane Nextseq sequencing ###
Overview
- A series of scripts that gets a consensus file, coverage file, and variant file for SARS-CoV-2 sequenced on Illumina nextseq with multiple libraries per lane
- copies files from raw data file on turbo to local directory based on the sample names in the sample sheet.

Requirements
- mkdir for each Library (or Sequencing Plate)
- add map file  with two colums (Processing_Plate sample_name)
	-in Sample_ACCN, replace _ with -
	-  map file shoul be named "Library".csv. Map will be used to move fastq files from raw_data folder.
- change "Library Name" to the name of the plate/library, and change "sequencing_run_path" to the raw data path for the sequencing run  in SARS-CoV-2_Illumina_nextseq.sh

Usage: sbatch SARS-CoV-2_Illumina_nextseq.sh

