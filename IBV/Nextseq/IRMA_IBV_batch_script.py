
import	argparse
from glob import glob
import pandas as pd
import subprocess
import os

parser= argparse.ArgumentParser ()
parser.add_argument('--prefix', action="store", dest='prefix')
args = parser.parse_args()

scripts_dir = "/nfs/turbo/umms-alauring/shared_projects/IBV_illumina/"+ args.prefix + "/scripts"
working_dir = "/nfs/turbo/umms-alauring/shared_projects/IBV_illumina/"+ args.prefix
trim_dir = "/nfs/turbo/umms-alauring/shared_projects/IBV_illumina/"+ args.prefix + "/data/cutadapt"
fastq_dir = "/nfs/turbo/umms-alauring/shared_projects/IBV_illumina/"+ args.prefix + "/data/fastq_renamed"

indivs = []
for filename in glob ("data/fastq_renamed/*.1.fastq.gz"):
	filename_only = filename.split("/")[2] 
	name = filename_only.split(".")[0]
	if not os.path.exists(name +"/*.fasta"):
		indivs.append(name)
print(indivs)


### bulk of sbatch script
def write_batch_script(indiv):
	sbatch_script = '''#!/bin/bash
#
#SBATCH --job-name=IBV_Illumina
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20000m 
#SBATCH --time=3:10:00
#SBATCH --account=alauring99
#SBATCH --partition=standard


cd {working_dir}

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -q 25 -m 20 -o  {trim_dir}/{indiv}.1.trimmed.fastq  \
-p  {trim_dir}/{indiv}.2.trimmed.fastq  {fastq_dir}/{indiv}.1.fastq.gz  {fastq_dir}/{indiv}.2.fastq.gz

IRMA FLU {trim_dir}/{indiv}.1.trimmed.fastq {trim_dir}/{indiv}.2.trimmed.fastq {indiv}

'''.format(trim_dir=trim_dir,indiv=indiv, fastq_dir=fastq_dir, working_dir=working_dir)
	return sbatch_script

## write sbatch script    
for i in indivs:
	sbatch_file = os.path.join(scripts_dir,i+'_irma.sbatch')
	indiv=i
	sbatch_script = write_batch_script(i)
    
	with open(sbatch_file, "w") as text_file:
		text_file.write(sbatch_script)

## submit sbatch script
sbatch_files = glob(os.path.join(scripts_dir,'*irma.sbatch'))
print(str(len(sbatch_files)) + ' files to submit.')

sbatch_files = glob(os.path.join(scripts_dir,'*irma.sbatch'))
print(str(len(sbatch_files)) + ' files to submit.')
for sf in sbatch_files:
	subprocess.call(['sbatch', sf])
