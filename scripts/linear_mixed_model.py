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

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', module='statsmodels')



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

    df = cf_analysis_lib.compatible_columns(df)
    metadata = cf_analysis_lib.compatible_columns(metadata)

    encoded_metadata = cf_analysis_lib.categories_to_numeric(metadata)
    encoded_metadata = cf_analysis_lib.remove_highly_correlated_data(encoded_metadata, 0.9, True)
    df = cf_analysis_lib.remove_highly_correlated_data(df, 0.9999, verbose=True)

    return df, encoded_metadata


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
        try:
            to_drop = list(df_combined_na.loc[:,df_combined_na.sum(axis=0) <1].columns)
        except TypeError as e:
            print(f"Error: {e} in iteration {i}", file=sys.stderr)
            for c in df_combined_na.columns:
                print(f"{c}: {df_combined_na[c].dtype}", file=sys.stderr)
            sys.exit(1)
        # fourth, drop any columns without enough unique values
        # we can't do this because it eliminates all the CS_ columns that are only 0 or 1
        # rcn = df_combined_na[subset_predictors].nunique()
        # to_drop += list(rcn[rcn < 10].index)

        if len(to_drop) > 0:
            df_combined_na = df_combined_na.drop(columns=to_drop)
        
        # fourth, make sure all the predictors are still present
        updated_predictors = list(set(subset_predictors).intersection(df_combined_na.columns))
        
        if include_culture_states:
            updated_culture_states = culture_states.intersection(set(df_combined_na.columns))
            formula = f"{dependent} ~ {' + '.join(updated_predictors)} + {' + '.join(updated_culture_states - {dependent})}"
        else:
            formula = f"{dependent} ~ {' + '.join(updated_predictors)} "

        result = None
        try:
            model = smf.mixedlm(
                formula=formula,
                data=df_combined_na,
                groups=df_combined_na["pwCF_ID"],
                method='BFGS'
            )
            result = model.fit()

        except Exception as e:
            print(f"Iteration {i} has error {e}\nformula: {formula}", file=sys.stderr)
            if isinstance(e, NameError):
                print(" ".join(list(df_combined_na.columns)))
                sys.exit(1)
            if isinstance(e, np.linalg.LinAlgError):
                print(f"LinAlgError {e} in iteration {i}", file=sys.stderr)
                for opt_method in 'Nelder-Mead', 'CG', 'COBYLA', 'Powell', 'trust-constr', 'trust-ncg':
                    try:
                        model = smf.mixedlm(
                            formula=formula,
                            data=df_combined_na,
                            groups=df_combined_na["pwCF_ID"],
                            method='BFGS'
                        )
                        result = model.fit(method=opt_method)
                        print(f"Success with {opt_method}", file=sys.stderr)
                        break
                    except Exception as e:
                        print(f"Error with {opt_method}: {e}", file=sys.stderr)

        if result is None:
            continue
        df_result = pd.DataFrame({
            "Predictor": result.params.index,
            "Estimate": result.params.values,
            "Std Err": result.bse.values,
            "P-value": result.pvalues.values
        })
        df_result["Iteration"] = i  # Add iteration number for tracking

        # Append to results
        results.append(df_result)

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

    all_predictors = df.columns
    dependent = re.sub(r'\W+', '_', args.dependent)

    if dependent not in encoded_metadata.columns:
        print(
            f"Error: {dependent} not in metadata so we can't predict it (check if it was dropped with --verbose)",
            file=sys.stderr)
        sys.exit(1)

    df_combined = df.merge(encoded_metadata, left_index=True, right_index=True, how='inner')
    df_combined = df_combined.dropna(subset=[dependent, "pwCF_ID"])

    if args.verbose:
        print(f"Running lmm with {len(all_predictors)} predictors", file=sys.stderr)
        print(f"Dependent variable is {dependent} with type {df_combined[dependent].dtypes}", file=sys.stderr)


    lmm(df=df_combined, dependent=dependent, all_predictors=all_predictors, num_predictors_per_model=args.predictors,
        num_iterations=args.iterations, output_file=args.output_file, include_culture_states=args.culture_states, verbose=args.verbose)
