### Author: Emily Bendall

### Purpose: concatenate individual consensus files and filters for completeness

# ================================= Run options ================================


# Usage:


"""

python ../filter_consensus.py \

    --prefix 20210712_Nanopore_Run_33 \
    
    --mode test \
   
    --min 13655

"""

# ============================= Import modules =================================

import	argparse
import	glob
import pandas as pd
import subprocess
from Bio import SeqIO
import numpy as np
import os


from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord

allowed_bases = ["A", "T", "C", "G", "a", "t", "c", "g"]

# ================================ Main =========================================

def main():

	parser= argparse.ArgumentParser ()
	parser.add_argument('--prefix', action="store", dest="prefix")
	parser.add_argument('--mode', action = "store", dest= "mode")
	parser.add_argument ('--min-length', action="store", dest="min", type=int)

	args = parser.parse_args()
# concat consesus sequences 

	collect_cmd = "cat */amended_consensus/*.fa > " + args.prefix + ".full.consensus.fa"

	if args.mode == "run":

		subprocess.call(collect_cmd, shell = True)

# filter for 90% completion
    
	if args.mode == "run":
		all_fasta = list ()
		full_consensus_file =  args.prefix +".full.consensus.fa"
		tmp_filename = full_consensus_file + ".tmp"
		final_filename =  args.prefix + ".90.consensus.fasta"

		for record in SeqIO.parse(full_consensus_file, "fasta"):

			sequence = str(record.seq)
			sequence_length = len ([b for b in sequence if b in allowed_bases])

			if (sequence_length >= args.min):
				all_fasta.append(record)
			else:
				print ("Sequence" + record.id + "was too short! Excluding from final file.")
             
			with open(tmp_filename, 'w') as full_fasta:
				SeqIO.write (all_fasta, full_fasta, "fasta")
      
		subprocess.call("sed '/^>/ s/ .*//' " + tmp_filename + " > " + final_filename, shell = True)
		subprocess.call("rm " + tmp_filename, shell = True)
		
# rename coverage files, add sample name, and concatenate

	cov_files =  "*/tables/*-coverage.txt"
	cov_filename = args.prefix+ ".coverage.csv"
	df_list = []
	
	for filename in glob.glob (cov_files):
	# change coverage file name per individual and add sample name
		path = filename.split("/")
		sample_name = path [0]
		df = pd.read_table(filename)
		df.columns =["chr", "pos", "cov", "consensus", "deletions", "ambiguous", "consensus_count", "consensus_average_quality"]
		df['ID'] = sample_name
		new_name =  sample_name + "/"+ sample_name + ".coverage.csv"
		df.to_csv(new_name)
		df_list.append(df)
	
	#concatenate coverage files
	df_full = pd.concat(df_list, axis=0, ignore_index=True)
	if args.mode == "run":
		df_full.to_csv(cov_filename, index = False)

if __name__ == "__main__":
	main ()
