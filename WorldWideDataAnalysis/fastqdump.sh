#!/bin/bash

#PBS -l ncpus=6
#PBS -l mem=8GB
#PBS -l walltime=10:00:00  
#PBS -N fasterqdump
#PBS -l jobfs=5GB
#PBS -P ob80
#PBS -l wd
#PBS -l other=hyperthread  
#PBS -l storage=scratch/ob80+gdata/ob80 



eval "$(conda shell.bash hook)"
conda activate sra-tools

mkdir -p $PBS_JOBFS/fasta fasta

for SRA in $(cat cf_metagenomes.txt); do 
	
	fasterq-dump $SRA --outdir $PBS_JOBFS/fasta --fasta-unsorted
	mv $PBS_JOBFS/fasta/* fasta/
done


