#!/bin/bash
  
#SBATCH --job-name=RSVA_ONT
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2000 
#SBATCH --time=3:10:00
#SBATCH --account=alauring99
#SBATCH --partition=standard

module load Bioinformatics cutadapt minimap2 samtools htslib/1.9-pcw6wxv ivar python/3.9.12 
mkdir data
mkdir data/fastq
python ../change_names_ONT.py --map $1 --mode run 
snakemake -s ../RSVA_ONT_Snakemake --cores 4 --rerun-incomplete --config run_name=$2
