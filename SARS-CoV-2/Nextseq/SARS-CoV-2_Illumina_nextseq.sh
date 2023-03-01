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
#python cp_fastq_files.py --sequencing_run "sequencing_run_path" --plate "Library Name"  --mode run
#python change_names_miseq.py -s "Library_name"/data/fastq -f "Library_name"/data/fastq_renamed
#module load python/3.9.1
#snakemake -s Snakemake-BWA_nextseq  --rerun-incomplete  --cores 4 --config run_name="Library Name"

## Script ##
module load Bioinformatics picard-tools R/4.2.0 ivar samtools fastqc  bwa bedtools2  htslib/1.9 
module load python/3.9.12
python cp_fastq_files.py --sequencing_run /nfs/turbo/umms-alauring/raw_data/2023/20230117_SC2_NS_R1  --plate 20220926_IAV_Illumina_Run_28 --mode run
python change_names_nextseq.py -s 20221219_SC2_Illumina_Run_73/data/fastq -f 20221219_SC2_Illumina_Run_73/data/fastq_renamed -run
#module load python/3.9.12 
snakemake -s Snakemake-BWA  --rerun-incomplete  --cores 4 --config run_name=20221219_SC2_Illumina_Run_73
