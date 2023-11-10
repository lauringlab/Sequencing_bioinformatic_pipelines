

# ============================= Import modules =================================

import argparse
import glob
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import glob

# ================================ Main =========================================


def main():


    parser = argparse.ArgumentParser()

    parser.add_argument('--map', action="store", dest="map") # mapping barcode to sample name. Header should be "Barcode" and "Name", with barcodes listed as NBXX.

    parser.add_argument('--mode', action="store", dest="mode") # run or test
    args = parser.parse_args()
    
    if args.mode== "run":
        subprocess.call ("mkdir data", shell=True)
        subprocess.call ("mkdir data/fastq", shell=True)

    # Read in barcode map and get barcode names to be used

    map = pd.read_csv(args.map, index_col = None, header = 0, dtype = object)
    
    map.columns = map.columns.astype(str)
    print (map.columns)

    barcodes = map.barcode
    sampleid = map.sample_id

    for bc in barcodes: 
        
        barcodes_list = barcodes.tolist()
        sampleid_list = sampleid.tolist()
        barcode_index = barcodes_list.index(bc)
        new_name = sampleid_list[barcode_index]

        if args.mode == "run":
             subprocess.call("cat fastq_pass/" + bc + "/*.fastq.gz > data/fastq/" +  new_name + ".fastq.gz", shell = True)

if __name__ == "__main__":

    main()

