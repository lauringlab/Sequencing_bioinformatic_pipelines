#import sys
import os
import argparse
import shutil
import gzip
import glob

parser = argparse.ArgumentParser(description='This program takes Nextseq fastq files and renames them as sample.read_direction.#.fastq and keeps a log of the change.')
parser.add_argument('-s',action='store',dest='s',help='The sorce directory containing the original fastq files')
parser.add_argument('-f',action='store',dest='f',help='The final directory that will hold the renamed fastq files')
parser.add_argument('-run',action='store_true',dest='test',default=False,help='Boolean switch to run program, without this the program runs in test mode: the log is made but no files are renamed')
#parser.add_argument('-gz',action='store_true',dest='zip',default=False,help='Boolean switch. Activate when working with .gz files')

args=parser.parse_args()
s=args.s
f=args.f
test=args.test
#zip=args.zip
# input argument is fastq directory

if not os.path.exists(f):
    os.makedirs(f)

# add slash to final location if it was forgotten
if f[-1] != '/':
    f=f+'/'
#print f

if s[-1] != '/':
    s=s+'/'
#print s

outfile = open(f+'renaming_log.txt','w')
#os.chdir(sys.argv[1])
if test==False:
    print ("running in test mode add option -run to run")
# copy fastq files
# The goal here is just to take out the junk added to the fastq file name
for zipfilename in glob.glob(s + "*.fastq.gz"):
    zipfilename= os.path.basename(zipfilename)
    name=zipfilename.split("_R")
    good_name = name[0].split("_")[0]
    read_number = name[1][0]

    perfect_name= good_name.replace("-","_")+"."+read_number+".fastq.gz"

    # Write file to new name

    print("COPYING "+ s+zipfilename + " to "+f+perfect_name)
    outfile.write(s+zipfilename + "\t COPIED to \t" + f+perfect_name + "\n")
    if test==True:
        shutil.copy(s+zipfilename,f+perfect_name)

outfile.close()

