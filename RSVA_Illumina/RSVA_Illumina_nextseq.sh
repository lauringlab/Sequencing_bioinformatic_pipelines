#!/bin/bash

#SBATCH --job-name=RSVA_Illumina
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4g 
#SBATCH --time=2:10:00
#SBATCH --account=alauring99
#SBATCH --partition=standard


## EXAMPLE ##

#conda activate snakemake
#module load Bioinformatics picard-tools R/3.6.1 ivar samtools fastqc  bwa bedtools2
#python cp_fastq_files.py --sequencing_run "sequencing_run_path" --plate "Library Name"  --mode run
#python change_names_miseq.py -s "Library_name"/data/fastq -f "Library_name"/data/fastq_renamed
#module load python/3.9.1
#snakemake -s Snakemake-BWA_nextseq  --rerun-incomplete  --cores 4 --config run_name="Library Name"

PROJECTSET=$1
PLATENAME=$2

## Script ##
module load Bioinformatics picard-tools R/4.2.0 ivar samtools fastqc  bwa bedtools2  htslib/1.9 
module load python/3.9.12
python /nfs/turbo/umms-alauring/shared_projects/RSVA_Illumina/cp_fastq_files.py --sequencing_run /nfs/turbo/umms-alauring/raw_data/2023/$PROJECTSET  --plate $PLATENAME --mode run
python /nfs/turbo/umms-alauring/shared_projects/RSVA_Illumina/change_names_nextseq.py -s $PLATENAME/data/fastq -f $PLATENAME/data/fastq_renamed -run
#module load python/3.9.12 
snakemake -s Snakemake-RSVA-Cutadapt  --rerun-incomplete  --cores 4 --config run_name=$PLATENAME
