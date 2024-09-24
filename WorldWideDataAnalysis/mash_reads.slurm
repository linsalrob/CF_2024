#!/bin/bash
#SBATCH --account=pawsey1018
#SBATCH --job-name=mash
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -o mash-%j.out
#SBATCH -e mash-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"

# Create a temporary environment to house mash and install it.
TMP=$(for i in {1..12}; do printf "%x" $((RANDOM % 16)); done); TMP="mash_$TMP";

echo "Creating mash in /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP"
mamba create -y --prefix=/scratch/pawsey1018/edwa0468/software/miniconda3/$TMP python=3.12 mash
conda activate /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP

FASTA="fasta"
MASH="mash"

for F in $(find $FASTA  -type f -printf "%f\n" | perl -ne 'm/^(...\d+)/; print "$1\n"' | awk '!s[$0]++'); do 
	if [[ ! -e $MASH/$F.msh ]]; then
		echo $F; 
		mash sketch -p 8 -I $F -r -o $MASH/$F $FASTA/${F}*;
	fi
done

rm -rf /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
