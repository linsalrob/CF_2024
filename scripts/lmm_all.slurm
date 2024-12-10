#!/bin/bash
#SBATCH --job-name=LinearMixedModelling
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH -o lmm-%A_%a.out
#SBATCH -e lmm-%A_%a.err

eval "$(conda shell.bash hook)"
conda activate bioinformatics
export PYTHONPATH=$PYTHONPATH:$HOME/GitHubs/CF_Data_Analysis

VAR="IP_vs_OP Hospital Age Gender Sample_Type H2_Uncorrected CH4_Uncorrected CO2 CH4_H2_ratio_corrected Corr_ NTM IgE Spec_IgE Precipitins FVC FEV1_best_FEV1 FEV1_RATIO_SCORE FEV1_Obstruction_Rank Cystic_Fibrosis_related_diabetes_CFRD_ Pancreatic_insufficiency_PI_ CF_gene_1 CFLD CS_mucoid CS_non_mucoid CS_Pseudomonas_aeruginosa CS_Oral_flora CS_Stenophotomonas_maltophilia CS_Aspergillus_fumigatus CS_Aspergillus_flavus CS_Candida_albicans CS_Mycobacteroides_abscessus CS_Mycobacterium_intracellulare CS_Staphylococcus_aureus CS_Inquilinus_limosus CS_Achromobacter_xylosoxidans CS_Burkholderia_cepacia CS_NTM__Smear_negative_ CS_Mycolicibacter_terrae CS_Aspergillus_nidulans CS_MAC__Smear_negative_ CS_Penicillium CS_Aspergillus_niger CS_Scedosporium_apiospermum CS_Lomentospora_prolificans CS_Acremonium_species CS_MDR_Pseudomonas_aeruginosa CS_Haemophilus_influenzae CS_Scedosporium_apiospermum_1 P12M_0_Nothing P12M_1_Staphylococcus_aureus P12M_2_Pseudomonas_aeruginosa P12M_3_Stenotrophomonas_maltophilia P12M_4_Rapid_NTM P12M_5_Slow_NTM P12M_6_Aspergillus P12M_7_Burkholderia_cepacia P12M_8_Haemophilus_influenza P12M_9_Streptococcus_pneumoniae P12M_10_Achromobacter_xylosoxidans P12M_11_Serratia_marcesens P12M_12_Inquilinus_limosus P12M_13_Scedosporum_apiospermum P12M_14_Moraxella_catarrhalis Positive_for_PA_within_next_12_months N12M_Pseudomonas_aeruginosa N12M_Pseudomonas_aeruginosa_multi_resistant_ N12M_mucoid N12M_non_mucoid N12M_MAC N12M_M_intracellulare N12M_Candida_albicans N12M_Aspergillus_fumigatus N12M_Stenotrophomonas_maltophilia N12M_Staphylococcus_aureus N12M_M_gordonae N12M_M_chelonae O_Penicillium_species O_Enterobacter_cloacae O_Aspergillus_Flavus O_Mycobacterium_abscessus O_Cupriavidus_gilardii O_Mycobacterium_intracellular O_Scedosporium_apiospermum O_Pseudomonas_aeruginosa_multi_resistant_ O_Streptococcus_pyogenes O_Pseudomonas_aeruginosa_non_mucoid_ O_Chryseobacterium_indologenes O_Bordetella_species O_Klebsiella_oxytoca O_Mycobacterium_mucogenicum O_Mycobacterium_scotochromogen O_Acinetobacter_species O_Enterobacter_species O_Aspergillus_niger O_Aspergillus_terreus O_Pseudomonas_aeruginosa_mucoid_ O_Proteus_mirabilis O_talaromyces_purpurogenus O_Trichosporon_mycotoxinivorans Antibiotics_YN IV INH PO two_Amphotericin_B_Ambisome__INH four_Amikacin_IV two_Amikacin_INH two_Amoxicillin_Potassium_clavulanate_Aug_Duo__PO three_Azithromycin_IV two_Azithromycin_PO three_Aztreonam_IV three_Cefopime_IV four_Cefoxitin_IV three_Ceftazidime_IV two_Ceftazidime_INH one_Cephalexin_PO two_Ciprofloxacin_PO two_Clarithromycin_PO two_Clofazimine_PO four_Colistin_IV two_Colistin_IHN one_Flucloaxcillin_PO three_Imipenem one_Itraconazole_Lozenoc__PO three_Ivacaftor_Kalydeco_ three_Meropenem_IV three_Methylpredinosolone_IV three_Omalizumab_SC three_piperacillin_sodium_tazobactam_sodium_Tazocin__IV two_prednisolone_PO one_Sulfamethoxazole_trimethoprim_Bactrim__PO three_tobramycin_IV two_tobramycin_INH Sum_of_antibiotics Sum_of_antifungals Sum_of_steroid_mabs DNA_extraction__conc DNA_Conc_ng_ul_ Mean_Size_BP Total_Clusters_Passing_Filter_Million_"
# create an array of variables
VAR_ARR=($VAR)
# notes:
# to print the elements:   echo ${VAR_ARR[@]}
# to print the number of elements:  echo "Number of elements: ${#VAR_ARR[@]}"
# this list has 152 elements

# choose one variable based on our task id
DEP=${NAMES_ARRAY[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "$DEP started at " `date`;
DSTART=`date +%s`;
python /home/edwa0468/GitHubs/CF_Data_Analysis/scripts/linear_mixed_model.py -d $DEP -o $DEP.tsv -p 100 -i 5000 --verbose --datadir '../../..' --culture_states
DEND=`date +%s`;
echo "$DEP took" $((DEND-DSTART)) "seconds";
