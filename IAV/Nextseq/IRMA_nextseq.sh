

## EXAMPLE
	#module load python/3.9.12
	#module load Bioinformatics 
	#module load cutadapt
	#module load R/4.2.0
	#cd "Plate name"
	#python ../cp_fastq_files.py --sequencing_run "Sequencing run" --map "Plate map" --mode run
	#python ../change_names_nextseq.py -s data/fastq -f data/fastq_renamed -run 
	#python ../IRMA_IAV_batch_script.py --prefix "Plate name" --mode run
# to submit save file and then command bash IRMA_nextseq.sh



module load python/3.9.12
module load Bioinformatics cutadapt R/4.2.0
cd 20221209_IAV_Illumina_Run_31
python ../cp_fastq_files.py --sequencing_run 20230124_IAV_NS_Run_1 --map 20221209_IAV_Illumina_Run_31.csv --mode run
python ../change_names_nextseq.py -s data/fastq -f data/fastq_renamed -run 
python ../IRMA_IAV_batch_script.py --prefix 20221209_IAV_Illumina_Run_31 --mode run
	
	
	