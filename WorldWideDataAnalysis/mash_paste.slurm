#!/bin/bash
#SBATCH --account=pawsey1018
#SBATCH --job-name=mashpaste
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH -o mash_paste-%j.out
#SBATCH -e mash_paste-%j.err

# set -euo pipefail
eval "$(conda shell.bash hook)"

############ README
#
# We can create a tmp directory, but with pawsey you will run out of space 
# as an alternative you can create the tmp directory and provide is on the sbatch command
# 

if [[ $# == 0 ]]; then
	echo "Making tmp mash installation. You might want to reconsider this if you run several jobs"

	# Create a temporary environment to house mash and install it.
	TMP=$(for i in {1..12}; do printf "%x" $((RANDOM % 16)); done); TMP="mash_$TMP";

	echo "Creating mash in /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP"
	mamba create -y --prefix=/scratch/pawsey1018/edwa0468/software/miniconda3/$TMP python=3.12 mash
	conda activate /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP

else
	conda activate $1
fi

OUTPUT="worldwide"
MASH="mash"

mash paste $OUTPUT $(find $MASH -type f)


if [[ $# == 0 ]]; then
	rm -rf /scratch/pawsey1018/edwa0468/software/miniconda3/$TMP
fi

