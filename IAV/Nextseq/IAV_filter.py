### Purpose: Get consensus genomes from Illumina sequence data for Influenza A.


# ================================= Run options ================================


# Usage:


"""

python ../Run_IRMA_Pipeline_IAV_illumina.py \

    --prefix 20210712_Nanopore_Run_33 \
    
    --mode test \
   
    --min 11819

"""

# ============================= Import modules =================================

import	argparse
import	glob
import pandas as pd
import subprocess
from Bio import SeqIO

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
	
	IDS = []
	for filename in glob.glob (args.prefix + "/data/fastq_renamed/*.1.fastq.gz"):
		filename_only = filename.split("/")[3] 
		name = filename_only.split(".")[0]
		IDS.append(name)
	for id in IDS:
		print (id)
		
	## Collect, rename, and filter by completeness
	
	print ("\n\n### (3/3) Organize, Filter, and Export ###\n\n")
	
	# collect together in one file
	
	for id in IDS:
		collect_cmd = "cat " + id + "/amended_consensus/*.fa > " + id + "/" + id + ".consensus.fasta"
		
		print (collect_cmd)	
		if args.mode == "run":
			subprocess.call(collect_cmd, shell=True)
			
	#create a single consensus sequence for each sample
	
	for id in IDS:
		collect_cmd1 = "grep -v '>' " +args.prefix + "/"+ id + "/"+ id + ".consensus.fasta >" + args.prefix + "/"+ id + "/"+ id + ".consensus.tmp" 
		collect_cmd2= "tr -d '\\n' < " +args.prefix + "/"+ id + "/"+ id + ".consensus.tmp > " +args.prefix + "/" + id + "/"+ id + ".consensus.tmp2"
		collect_cmd3= "sed  -i '1i >" +id + "' " +args.prefix + "/"+ id + "/"+ id + ".consensus.tmp2" 
		collect_cmd4= "sed -i -e '$a\\' " +args.prefix + "/" + id + "/"+ id + ".consensus.tmp2"

		if args.mode == "run":

			subprocess.call(collect_cmd1, shell = True)
			subprocess.call(collect_cmd2, shell = True)
			subprocess.call(collect_cmd3, shell = True)
			subprocess.call(collect_cmd4, shell = True)
			subprocess.call("rm " +args.prefix + "/"+ id + "/"+ id + ".consensus.tmp" , shell = True)
            
	# concat consesus sequences 

	collect_cmd = "cat " + args.prefix + "/*/*.consensus.tmp2 > " +args.prefix + "/" + args.prefix + ".full.consensus.fa"

	if args.mode == "run":

		subprocess.call(collect_cmd, shell = True)
         
	# filter for 90% completion
    
	if args.mode == "run":
		all_fasta = list ()
		full_consensus_file = args.prefix + "/"+ args.prefix + ".full.consensus.fa"
		tmp_filename = full_consensus_file + ".tmp"
		final_filename = args.prefix + "/"+ args.prefix + ".90.consensus.fasta"

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
         
         
	# rename individual sample segments sequence
	for id in IDS:        
		for filename in glob.glob(args.prefix + "/"+ id + "/A_*.fasta"):
        
			if args.mode == "run":
				all_fasta = list ()
				tmp_filename = filename + ".tmp"
                
      
				for record in SeqIO.parse(filename, "fasta"):
					segment = str(record.id).split("_")[1]
					new_ID = id + "_" + segment
					record.id = new_ID
					all_fasta.append(record) 
  
          
				name = filename.split("/")[2]
				final_filename = args.prefix + "/"+ id + "/" + id  + "_" +name 
                

				with open(tmp_filename, 'w') as full_fasta:
					SeqIO.write (all_fasta, full_fasta, "fasta")


				subprocess.call("sed '/^>/ s/ .*//' " + tmp_filename + " > " + final_filename, shell = True)
				subprocess.call("rm " + tmp_filename, shell = True)
                
    #collect HA segments into subtype folders 

	collect_cmd_H3 = "cat " +args.prefix +  "/*/*_A_HA_H3.fasta > "  +args.prefix + "/"+ args.prefix + ".HA.H3.consensus.fa"
	collect_cmd_H1 = "cat " +args.prefix + "/*/*_A_HA_H1.fasta > " +args.prefix + "/" + args.prefix + ".HA.H1.consensus.fa"

	if args.mode == "run":

		subprocess.call(collect_cmd_H3, shell = True)
		subprocess.call(collect_cmd_H1, shell = True)

	if args.mode == "run":
		for filename in  glob.glob(args.prefix + "/" + args.prefix +".HA.H*.consensus.fa"):
			all_fasta = list ()
			tmp_filename = filename + ".tmp"
			final_filename = filename 
                
      
			for record in SeqIO.parse(filename, "fasta"):
				fasta_record = str(record.id).split("_")[0]
				record.id = fasta_record
				all_fasta.append(record)  

			with open(tmp_filename, 'w') as full_fasta:
				SeqIO.write (all_fasta, full_fasta, "fasta")


			subprocess.call("sed '/^>/ s/ .*//' " + tmp_filename + " > " + final_filename, shell = True)
			subprocess.call("rm " + tmp_filename, shell = True)
    
    # make new directory and copy individual segment sequences to the folder
	collect_cmd = "mkdir " +args.prefix + "/Segment_sequences"
	if args.mode == "run":

		subprocess.call(collect_cmd, shell = True)
	
	for id in IDS:
		collect_cmd = "cp " + args.prefix + "/"+ id + "/" + id + "_*.fasta " +args.prefix + "/Segment_sequences"
		print (collect_cmd)

		if args.mode == "run":

			subprocess.call(collect_cmd, shell = True)
            
	print("\n\n### Congratulations, the pipeline is complete! ###\n\n")
   # remove empty files
	remove_cmd = "find . -maxdepth 1 -type f -empty -print -delete"
	if args.mode == "run":
		subprocess.call (remove_cmd, shell=True)	
if __name__ == "__main__":

    main()