# ============================= How to run this pipeline ==========================


#  Create barcode map run_name_
	# two columns (Processing_Plate	sample_name)
	# in Sample_ACCN, replace _ with -
	

#  Run pipeline with usage as below.


# ================================= Run options ================================


# Usage:


# ============================= Import modules =================================

import argparse
import os
import subprocess
import pandas as pd

# ================================ Main =========================================

def main():

	parser= argparse.ArgumentParser ()
	parser.add_argument('--sequencing_run', action="store", dest="sequencing")
	parser.add_argument('--mode', action = "store", dest= "mode")
	parser.add_argument('--map', action="store", dest="map") # mapping barcode to sample name. Header should be "Barcode" and "Name", with barcodes listed as NBXX.
	args = parser.parse_args()
	# Read in barcode map and get barcode names to be used

	map = pd.read_csv(args.map, index_col = None, header = 0, dtype = object)
	map.columns = map.columns.astype(str)
	barcodes = map.sample_name
	
	mkdir_cmd = "mkdir data \n mkdir data/fastq \n mkdir data/cutadapt"
	if args.mode == "run":
		subprocess.call (mkdir_cmd, shell = True)
		
	for bc in barcodes:
		cp_cmd = "cp /nfs/turbo/umms-alauring/raw_data/2023/"  + args.sequencing + "/" + bc + "*.fastq.gz data/fastq"
		
		if args.mode == "run":
			subprocess.call (cp_cmd, shell = True)

if __name__ == "__main__":

    main()

	
