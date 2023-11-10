# Author: Andrew Valesano
# Purpose: Combine the consensus files and coverage data.

# Usage: python ivar_CombineAndExport.py --run-info test --min-length 100

# ======================= Import modules ======================

import argparse
import glob
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--min-length', action="store", dest="min", type = int)
    parser.add_argument('--run-name', action="store", dest="run")
    args = parser.parse_args()

    ### Edit and combine the consensus files into one fasta with 90% completeness ###
    consensus_files =  "data/ivar_consensus/*.fa"
    all_fasta = list()
    all_fasta_full = list ()
    tmp_filename =  "data/" +  "tmp.consensus.fasta"
    tmp_filename2 ="data/" +  "tmp.consensus2.fasta"
    final_filename =  args.run + ".90.consensus.fasta"
    final_filename_all =  args.run + ".full.consensus.fasta"

    for file in glob.glob(consensus_files):

        print("Working on file: " + file)
        filename_only = file.split("/")[2]
        id = filename_only.split(".")[0]
        for record in SeqIO.parse(file, "fasta"):
            record.id = id
            length = sum(map(lambda x: record.seq.count(x), ["a", "t", "g", "c", "A", "T", "G", "C"]))
            if(length >= args.min):
                all_fasta.append(record)
            else:
                print("Sequence " + id + " was too short! Excluding from final file. Length was ")
        print ("gone through all files")
    with open(tmp_filename, 'w') as filtered_fasta:
        SeqIO.write(all_fasta, filtered_fasta, "fasta")
    os.system("sed '/^>/ s/ .*//' " + tmp_filename + " > " + final_filename)
    os.system("rm " + tmp_filename)

## combine all fasta files 

    for file in glob.glob(consensus_files):

        print("Working on file: " + file)
        filename_only = file.split("/")[2]
        id = filename_only.split(".")[0]
        for record in SeqIO.parse(file, "fasta"):
            record.id = id
            all_fasta_full.append(record)
    with open(tmp_filename2, 'w') as full_fasta:
        SeqIO.write(all_fasta_full, full_fasta, "fasta")
    os.system("sed '/^>/ s/ .*//' " + tmp_filename2 + " > " + final_filename_all)
    os.system("rm " + tmp_filename2)
### Combine coverage files ###
    cov_files = "data/coverage/*.csv"
    cov_filename = args.run + ".coverage.csv"

    df_list = []
    for filename in glob.glob(cov_files):
        if os.path.getsize(filename) > 0:
            filename_only = filename.split("/")[2]
            id = filename_only.split(".")[0]
            df= pd.read_csv(filename, sep = '\t')
            df.columns = ["chr", "pos", "cov"]
            df['ID'] = id
            df_list.append(df)
    df_full = pd.concat(df_list, axis=0, ignore_index=True)
    df_full['run'] = args.run
    df_full.to_csv(cov_filename, index = False)


if __name__ == "__main__":
    main()
