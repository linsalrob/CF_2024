#!/bin/bash
#SBATCH --job-name=fastqdump
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -o fqdump-%j.out
#SBATCH -e fqdump-%j.err

eval "$(conda shell.bash hook)"
conda activate sra

TMPFS=$BGFS

mkdir -p $TMPFS/fasta fasta

IFS=$'\n'
for SRALINE in $(cat cf_metagenome_locations.tsv); do 
	SRA=$(echo $SRALINE | awk '{print $1}');
	if [[ ! -e fasta/$SRA.fasta.gz ]]; then
		LOC=$(echo $SRALINE | awk '{print $2}');
		if [[ $LOC != "None" ]]; then
			fasterq-dump $SRA --outdir $TMPFS/fasta --fasta-unsorted -e 8
			pigz $TMPFS/fasta/$SRA.fasta
			mv $TMPFS/fasta/* fasta/
		fi
	fi
done



