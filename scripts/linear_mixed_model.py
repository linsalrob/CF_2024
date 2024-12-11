"""
A linear mixed model for the CF analysis pipeline.

This comes from the LMM.ipynb notebook, but is designed so you can run it as a script to allow lots of iterations
"""


__author__ = 'Rob Edwards'
__version__ = '0.1.0'

import os
import sys
import argparse
import re
import numpy as np
import pandas as pd
import random

import statsmodels.formula.api as smf
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

datadir = '..'
sys.path.append('..')
import cf_analysis_lib


def read_data_frames(sequence_type = "MGI", datadir = "..", sslevel = "subsystems_norm_ss.tsv.gz",  taxa = "genus", verbose=False):

    if verbose:
        print('Reading data frames', file=sys.stderr)

    ss_df = cf_analysis_lib.read_subsystems(
        os.path.join(datadir, sequence_type, "FunctionalAnalysis", "subsystems", sslevel), sequence_type)
    ss_df = ss_df.T
    if verbose:
        print(f"The subsystems df has shape: {ss_df.shape}", file=sys.stderr)

    genus_otu = cf_analysis_lib.read_taxonomy(datadir, sequence_type, taxa)
    genus_otu = genus_otu.T
    if verbose:
        print(f"The taxonomy df has shape: {genus_otu.shape}", file=sys.stderr)
    metadata = cf_analysis_lib.read_metadata(datadir, sequence_type, categorise=True)
    if verbose:
        print(f"The metadata df has shape: {metadata.shape}", file=sys.stderr)

    df = ss_df.merge(genus_otu, left_index=True, right_index=True, how='inner')

    # rename the columns to remove any special characters
    df.columns = [re.sub(r'\W+', '_', col) for col in df.columns]
    metadata.columns = [re.sub(r'\W+', '_', col) for col in metadata.columns]

    # rename the columns that start with numbers
    nos = {"1": "one", "2": "two", "3": "three", "4": "four", "5": "five", "6": "six", "7": "seven", "8": "eight",
           "9": "nine"}
    new_df_cols = {}
    for c in df.columns:
        if c[0] in nos:
            new_df_cols[c] = c.replace(c[0], nos[c[0]], 1)
    df = df.rename(columns=new_df_cols)

    newmetacols = {}
    for c in metadata.columns:
        if c[0] in nos:
            newmetacols[c] = c.replace(c[0], nos[c[0]], 1)
    metadata = metadata.rename(columns=newmetacols)

    # convert categories to numbers
    # we make a copy of each data frame and then overwrite the columns. Internally this is more efficient than making a df one
    # column at a time.
    encoded_metadata = metadata.copy()
    to_delete = ['Pseudomonas_Culture']  # this is a duplicate columns with 'CS_Pseudomonas_aeruginosa'
    for col in metadata:
        if pd.api.types.is_numeric_dtype(metadata[col]):
            continue
        elif col == sequence_type:
            continue
        elif metadata[col].dtypes == 'category':
            encoded_metadata[col] = metadata[col].cat.codes
        else:
            if verbose:
                print(f'Dropping {col}: {metadata[col].dtypes}', file=sys.stderr)
            to_delete.append(col)

    encoded_metadata = encoded_metadata.drop(columns=to_delete)

    return df, encoded_metadata

def remove_highly_correlated_variables(df, cutoff=0.9, verbose=False):
    if verbose:
        print("Calculating correlation matrix", file=sys.stderr)
    correlation_matrix = df.corr(method='pearson')

    # Find highly correlated pairs (absolute correlation > 0.8)
    high_corr_abundance = correlation_matrix.unstack().reset_index()
    high_corr_abundance.columns = ['From', 'To', 'Correlation']
    high_corr_abundance = high_corr_abundance[
        (high_corr_abundance['Correlation'].abs() > 0.9999) & (high_corr_abundance['From'] != high_corr_abundance['To'])
        ]

    # Drop duplicate pairs (e.g., (A, B) and (B, A))
    high_corr_abundance = high_corr_abundance.drop_duplicates(subset=['Correlation'])

    if verbose:
        print(f"Dropping {len(high_corr_abundance['To'])} highly correlated variables", file=sys.stderr)
    df = df.drop(columns=high_corr_abundance['To'])
    return df

def lmm(df, dependent, all_predictors, num_predictors_per_model=100, num_iterations=10000, output_file=None, include_culture_states=False, verbose=False):
    """
    Run a linear mixed model with a random subset of predictors
    :param df: The combined data frame
    :param dependent: the dependent variable
    :param all_predictors: the columns we should use as predictors (probably the OTUs and functional categories)
    :param num_predictors_per_model: how many predictors to use in each model
    :param num_iterations: how many models to run
    :return:
    """
    results = []

    culture_states = {'CS_mucoid', 'CS_non_mucoid', 'CS_Pseudomonas_aeruginosa', 'CS_Oral_flora',
                      'CS_Stenophotomonas_maltophilia', 'CS_Aspergillus_fumigatus', 'CS_Aspergillus_flavus',
                      'CS_Candida_albicans', 'CS_Mycobacteroides_abscessus', 'CS_Mycobacterium_intracellulare',
                      'CS_Staphylococcus_aureus', 'CS_Inquilinus_limosus', 'CS_Achromobacter_xylosoxidans',
                      'CS_Burkholderia_cepacia', 'CS_NTM__Smear_negative_', 'CS_Mycolicibacter_terrae',
                      'CS_Aspergillus_nidulans', 'CS_MAC__Smear_negative_', 'CS_Penicillium', 'CS_Aspergillus_niger',
                      'CS_Lomentospora_prolificans', 'CS_Acremonium_species',
                      'CS_MDR_Pseudomonas_aeruginosa', 'CS_Haemophilus_influenzae'}

    not_currently_in_metadata = set()
    for c in culture_states:
        if c not in encoded_metadata.columns:
            not_currently_in_metadata.add(c)
    culture_states = culture_states - not_currently_in_metadata

    for i in range(num_iterations):
        if verbose and i % (num_iterations/10) == 0:
            print(f"Iteration {i}", file=sys.stderr)
        
        # first, come up with a random list of predictors
        subset_predictors = random.sample(list(all_predictors), num_predictors_per_model)  # Randomly select predictors
        # second, drop any rows where the predictors have null values.
        # this will change the outcome of the next line!
        df_combined_na = df.dropna(subset=subset_predictors)
        # third, remove any columns whose column sum is 0
        to_drop = list(df_combined_na.loc[:,df_combined_na.sum(axis=0) <1].columns)
        # fourth, drop any columns without enough unique values
        rcn = df_combined_na[subset_predictors].nunique()
        to_drop += list(rcn[rcn < 10].index)

        if len(to_drop) > 0:
            df_combined_na = df_combined_na.drop(columns=to_drop)
        
        # fourth, make sure all the predictors are still present
        updated_predictors = list(set(subset_predictors).intersection(df_combined_na.columns))
        
        if include_culture_states:
            formula = f"{dependent} ~ {' + '.join(updated_predictors)} + {' + '.join(culture_states - {dependent})}"
        else:
            formula = f"{dependent} ~ {' + '.join(updated_predictors)} "

        try:
            with warnings.catch_warnings():
                # this just suppresses some statsmodels warnings.
                warnings.filterwarnings("ignore", module="statsmodels")
                model = smf.mixedlm(
                    formula=formula,
                    data=df_combined_na,
                    groups=df_combined_na["pwCF_ID"]
                )
                result = model.fit()

            df_result = pd.DataFrame({
                "Predictor": result.params.index,
                "Estimate": result.params.values,
                "Std Err": result.bse.values,
                "P-value": result.pvalues.values
            })
            df_result["Iteration"] = i  # Add iteration number for tracking

            # Append to results
            results.append(df_result)
        except Exception as e:
            print(f"Iteration {i} has error {e} formula: {formula}", file=sys.stderr)

    # Combine results into a single DataFrame
    all_results = pd.concat(results)

    # Group by predictor and aggregate results
    final_results = (
        all_results.groupby("Predictor")[["Estimate", "Std Err", "P-value"]]
        .mean()
        .rename(columns={"Estimate": "Mean Estimate", "Std Err": "Mean Std Err", "P-value": "Mean P-value"})
        .reset_index()
    )

    with open(output_file, 'w') as out:
        final_results.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the lmm lots of times')
    parser.add_argument('-o', '--output_file', help='ouput file', required=True)
    parser.add_argument('-d', '--dependent', help='dependent variable (note we will reformat the name)', required=True)
    parser.add_argument('-p', '--predictors', help='number of predictors per model', type=int, default=100)
    parser.add_argument('-s', '--sequence_type', help='sequence type', default='MGI')
    parser.add_argument('-l', '--sslevel', help='subsystem level', default='subsystems_norm_ss.tsv.gz')
    parser.add_argument('-n', '--datadir', help='data directory', default='..')
    parser.add_argument('-t', '--taxa', help='taxonomic level', default='family')
    parser.add_argument('-i', '--iterations', help='number of iterations', type=int, default=10000)
    parser.add_argument('-c', '--culture_states', help='include culture state (CS) metadata in the predictors', action='store_true')

    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    print(f"Running linear_mixed_model.py with arguments: {args}", file=sys.stderr)

    df, encoded_metadata = read_data_frames(sequence_type=args.sequence_type, sslevel=args.sslevel, taxa=args.taxa, datadir=args.datadir, verbose=args.verbose)

    df = remove_highly_correlated_variables(df, 0.9999)
    encoded_metadata = remove_highly_correlated_variables(encoded_metadata.drop(columns=args.sequence_type), 0.9, verbose=args.verbose)

    if args.dependent not in encoded_metadata.columns:
        print(
            f"Error: {args.dependent} not in metadata so we can't predict it (check if it was dropped with --verbose)",
            file=sys.stderr)
        sys.exit(1)

    df_combined = df.merge(encoded_metadata, left_index=True, right_index=True, how='inner')

    all_predictors = df.columns
    dependent = re.sub(r'\W+', '_', args.dependent)

    lmm(df=df_combined, dependent=dependent, all_predictors=all_predictors, num_predictors_per_model=args.predictors,
        num_iterations=args.iterations, output_file=args.output_file, include_culture_states=args.culture_states, verbose=args.verbose)
