#!/bin/bash

#SBATCH --job-name=SARS_CoV-2_Illumina
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=1500m 
#SBATCH --time=2:10:00
#SBATCH --account=alauring99
#SBATCH --partition=standard


## EXAMPLE ##

#conda activate snakemake
#module load Bioinformatics picard-tools R/3.6.1 ivar samtools fastqc  bwa bedtools2
#python change_names_miseq.py -s nfs/turbo/umms-alauring/raw_data/2022/Sequencing_Plate -f Sequencing_Plate/data/fastq_renamed
#module load python/3.9.1
#snakemake -s Snakemake-BWA  --rerun-incomplete  --cores 4 --config run_name=Sequencing_Plate

## Script ##
#conda activate snakemake
module load Bioinformatics picard-tools R/4.2.0 ivar samtools fastqc  bwa bedtools2  htslib/1.9 
#python change_names_miseq.py -s /nfs/turbo/umms-alauring/raw_data/2022/20221207_SC2_Illumina_Run_72 -f 20221207_SC2_Illumina_Run_72/data/fastq_renamed -run
module load python/3.9.12 
snakemake -s Snakemake-BWA  --rerun-incomplete  --cores 4 --config run_name=20221207_SC2_Illumina_Run_72
