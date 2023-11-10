README for RSVB Illumina Processing 

Overview
- For processing RSVB sequences derived from overlapping (~400bp) amplicons and sequenced on an Illumina Nextseq. 
- Primers are trimmed using primer sequences with cutadapt. Reads are aligned with BWA (gisaid reference EPI_ISL_1653999 -hRSV/B/Australia/VIC-RCH056/2019) and a consensus sequence is called using ivar.
- 2 python scripts to copy (cp_fastq_files.py) and rename (change_names_nextseq.py) necessary fastq files 
- Snakemake file that processes fastq files and creates consensus sequence


Requirements
- Need to have snakmake installed
- Need to have biopython and pandas installed in python/3.9.12

Steps
1. Make a folder for each sequencing plate inside RSVB_Illumina
2. Add barcode map to plate folder. Barcode map is a 2 column csv file with "Processing_Plate,sample_name" as the header
3. sbatch RSVB_Illumina_nextseq.sh "sequencing run name" "plate name"

 