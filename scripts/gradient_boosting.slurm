#!/bin/bash
#SBATCH --job-name=GradientBoosting
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH -o gb-%j.out
#SBATCH -e gb-%j.err

eval "$(conda shell.bash hook)"
conda activate bioinformatics
export PYTHONPATH=$PYTHONPATH:$HOME/GitHubs/CF_Data_Analysis


python /home/edwa0468/GitHubs/CF_Data_Analysis/scripts/gradient_boosting.py -r gb_1000_1000_family_ss.tsv -n 1000 -i gb_1000_1000_family_ss_images -e 1000 --sequence_type MGI --taxa family --verbose --datadir '../..' 
