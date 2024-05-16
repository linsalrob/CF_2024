#!/bin/bash
#PBS -P ob80
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l jobfs=1GB
#PBS -l walltime=60:00
#PBS -l wd
#PBS -l other=hyperthread
#PBS -l storage=scratch/ob80+gdata/ob80
#PBS -N fastp
#PBS -r y


###################### NOTES ########################
# 
# READ THIS: https://opus.nci.org.au/display/Help/nci-parallel
# NCI doesn't like you using array jobs!!
#
#
# You can only submit 9 jobs (1-10) per array, but you
# can easily set up a for loop to submit multiple array jobs
# Yes. I know, this defeats the point
#
# for R in $(seq 11 10 111); do  E=$((R+9)); echo "$R $E"; qsub -J $R-$E assemble_all.sh; done
#
#  NOTE: The submit command is -J 1-10 and you MUST include -r y in the commands
# NOTE2: The variable is called $PBS_ARRAY_INDEX


set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics

FN=$(head -n $PBS_ARRAY_INDEX R1_reads.txt | tail -n 1)
echo "PROCESSING $FN"

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

# this is a hack because one some clusters we have a FS called BGFS
# TMPFS=$BGFS
TMPFS=$PBS_JOBFS

cp /home/edwa0468/GitHubs/fast-adapter-trimming/adapters/IlluminaAdapters.fa $TMPFS

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}
cp $SOURCE/$R1 $SOURCE/$R2 $TMPFS/fastq

echo "fastp -n 1 -l 100 -i $TMPFS/fastq/$R1 -I $TMPFS/fastq/$R2 -o $TMPFS/output/$R1 -O $TMPFS/output/$R2 --adapter_fasta $TMPFS/IlluminaAdapters.fa " >&2
fastp -n 1 -l 100 -i $TMPFS/fastq/$R1 -I $TMPFS/fastq/$R2 -o $TMPFS/output/$R1 -O $TMPFS/output/$R2 --adapter_fasta $TMPFS/IlluminaAdapters.fa

mkdir -p fastq_fastp
mv $TMPFS/output/$R1 $TMPFS/output/$R2  fastq_fastp
