#!/bin/bash

#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=10:00:00  
#PBS -N mash
#PBS -l jobfs=5GB
#PBS -P ob80
#PBS -q normal
#PBS -l wd
#PBS -l other=hyperthread  
#PBS -l storage=scratch/ob80+gdata/ob80 



eval "$(conda shell.bash hook)"
conda activate mash

TMPFS=$PBS_JOBFS

mkdir -p $TMPFS/mash mash_sketches

R1=$(head -n $JOBNUM R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}
MSH=${R1/_R1.fasta.gz/}
echo "R1: $R1 R2: $R2 MASH: $MSH" >&2
if [[ -e cutadapt/$R2 ]]; then
	mash sketch -r -o $TMPFS/mash/$MSH cutadapt/$R1 cutadapt/$R2
else
	mash sketch -r -o $TMPFS/mash/$MSH cutadapt/$R1 
fi

# mash paste fasta.msh $TMPFS/mash/*
cp $TMPFS/mash/* mash_sketches
rm -rf $TMPFS/mash/

