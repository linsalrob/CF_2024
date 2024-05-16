#!/bin/bash
#SBATCH --job-name=splitfasta
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH -o splitfasta-%j.out
#SBATCH -e splitfasta-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics



split_fasta() {
	FA=$1;
	R1=${FA/.fasta.gz/_R1.fasta}
	R2=${FA/.fasta.gz/_R2.fasta}
	if [[ -e fasta_split/$R1.gz ]]; then return 0; fi
	/home/edwa0468/GitHubs/EdwardsLab/bin/fasta_split $BGFS/fasta/$FA $BGFS/fasta_split/$R1 $BGFS/fasta_split/$R2
	pigz $BGFS/fasta_split/$R1 
	pigz $BGFS/fasta_split/$R2
	mv -t fasta_split/ $BGFS/fasta_split/$R1.gz $BGFS/fasta_split/$R2.gz
}

export -f split_fasta

cp -r fasta $BGFS
mkdir -p $BGFS/fasta_split fasta_split

find -L $BGFS/fasta -type f -name \*fasta.gz -printf "%f\n" | parallel -j 16 split_fasta



