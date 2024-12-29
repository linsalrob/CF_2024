#!/bin/bash
#SBATCH --job-name=GradientBoosting
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH -o gb-%j.out
#SBATCH -e gb-%j.err

eval "$(conda shell.bash hook)"
conda activate bioinformatics
export PYTHONPATH=$PYTHONPATH:$HOME/GitHubs/CF_Data_Analysis

# usage: gradient_boosting.py [-h] -r RESULTS [-n ITERATIONS] [-e NESTIMATORS] [-f FEATURES] [-y IMAGE_FEATURES] [-i IMAGES] [-d DATADIR] [-s SEQUENCE_TYPE] [-l SSLEVEL] [-t TAXA] [-v]
# gradient_boosting.py: error: the following arguments are required: -r/--results

ITERATIONS=1000
ESTIMATORS=1000
TAXA=family
python /home/edwa0468/GitHubs/CF_Data_Analysis/scripts/gradient_boosting.py \
				-r gb_${ITERATIONS}_${ESTIMATORS}_${TAXA}_ss.tsv \
				-i gb_${ITERATIONS}_${ESTIMATORS}_${TAXA}_ss_images \
				-n $ITERATIONS \
				-e $ESTIMATORS \
				-y 20 \
				-f 100 \
				--sequence_type MGI \
				--taxa $TAXA \
				--verbose \
				--datadir '../..' \
				--printall
