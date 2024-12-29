"""
Read different data types into a pandas dataframe.

Includes support for:
 - subsystems
 - taxonomy
"""

import os
import sys
import pandas as pd
from sklearn.impute import SimpleImputer

from .metadata_data import metadata_definitions

corrections = {
    "MGI" : {
        '1112926_20171212_S': '1447437_20171212_S',
        '1128691_20170206_S': '1128691_20171206_S',
        '1255498_20171212_S': '1590009_20171212_S',
        '1316979_20171215_S': '1651490_20171215_S',
        '1598281_20180508_S': '1588281_20180508_S',
        '1723809_20180227_S': '1085876_20180227_S',
        '649354_20170206_S': '639354_20171206_S',
        '652927_20180226_S': '715927_20180226_S',
        '658355_20180301_S': '658355_20180327_S',
        '777851_20170918_S': '778851_20170918_S',
        '788707_20181126_S': '788707_20181129_S'
    },
    "minion" : {
        '1112926_20171212_S': '1447437_20171212_S',
        '1255498_20171212_S': '1590009_20171212_S',
        '1316979_20171215_S': '1651490_20171215_S',
        '1598281_20180508_S': '1588281_20180508_S',
        '698917_20190119_S': '698917_20180119_S'
        }
}


def read_taxonomy(datadir, sequence_type, taxonomy):
    """
    Read the taxonomy file and return a data frame
    """

    if sequence_type.lower() == 'mgi':
        sequence_type = 'MGI'
    elif sequence_type.lower() == 'minion':
        sequence_type = 'minion'
    else:
        raise ValueError(f"Sorry. Don't know what sequence type {sequence_type} is supposed to be")

    tax_file = os.path.join(datadir, sequence_type, "Taxonomy", f"{sequence_type}_reads_{taxonomy}.normalised.tsv.gz")
    if not os.path.exists(tax_file):
        raise FileNotFoundError(f"Error: {tax_file} does not exist")
    df = pd.read_csv(tax_file, sep='\t', compression='gzip')
    df = df[df['taxonomy'].str.contains('k__Bacteria')]
    df = df[~df['taxonomy'].str.endswith(f'{taxonomy[0]}__')]
    df = df.set_index('taxonomy')
    df = df.rename(columns=corrections[sequence_type])
    df.index = df.index.str.replace(f'{taxonomy[0]}__', '').str.replace('Candidatus ', '')
    df.index = df.index.str.split(';').str[-1]

    df = df.sort_index(axis=1)
    return df

def read_metadata(datadir, sequence_type, categorise=False, verbose=False):
    """
    Read the metadata file and return a data frame
    """
    sequencing = []
    if sequence_type.lower() == 'mgi':
        sequencing = ['MGI']
    elif sequence_type.lower() == 'minion':
        sequencing = ['minion']
    elif sequence_type.lower() == 'mgi_minion':
        sequencing = ['MGI', 'minion']
    else:
        raise ValueError(f"Sorry. Don't know what sequence type {sequence_type} is supposed to be")

    metadata_file = os.path.join(datadir, "Metadata", "Metadata.tsv")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Error: {metadata_file} does not exist")
    metadata = pd.read_csv(metadata_file, encoding='windows-1252', sep="\t", index_col=0)

    if len(sequence_type) == 1:
        metadata = metadata[~metadata[sequence_type[0]].isna()]

    metadata = metadata.rename(columns={'Pseudomonas': 'Pseudomonas Culture'})

    for ix in metadata.index:
        for seq_type in sequencing:
            s = metadata.loc[ix, seq_type]
            if s in corrections[seq_type]:
                metadata.loc[ix, seq_type] = corrections[seq_type][s]

    # impute missing values by most frequent (i.e. mode)
    # imputer = SimpleImputer(strategy='most_frequent')
    # mean_imputer = SimpleImputer(strategy='mean')
    if categorise:
        if verbose:
            print(f"Converting {c} to category. Make sure to dropna(axis=1)", file=sys.stderr)
        # convert the metadata to categories!
        mdx_types = metadata_definitions()
        for c in metadata.columns:
            if c in mdx_types and mdx_types[c] == 'Categorical':
                metadata[c] = metadata[c].astype('category')
            elif c in mdx_types and mdx_types[c] == 'Date':
                metadata[c] = pd.to_datetime(metadata[c])

    return metadata

def read_subsystems(subsystems_file, sequence_type):
    """
    Read the subsystems file and return a data frame
    """
    if not os.path.exists(subsystems_file):
        raise FileNotFoundError(f"Error: {subsystems_file} does not exist")

    if subsystems_file.endswith('.gz'):
        df = pd.read_csv(subsystems_file, sep='\t', compression='gzip', index_col=0)
    else:
        df = pd.read_csv(subsystems_file, sep='\t', index_col=0)
    df = df.rename(columns=corrections[sequence_type])
    return df

def sorted_presence_absence(df1, df2, minrowsum=0, asc_sort=False):
    """
    Join the two tables and return the sorted version
    """
    # filter so we only include samples sequenced on both MGI and MinION
    common_columns = df1.columns.intersection(df2.columns)
    df1_both = df1[common_columns]
    df2_both = df2[common_columns]

    # create a presence/absence matrix
    df1_presence = (df1_both > 0).astype(int)
    df2_presence = (df2_both > 0).astype(int)*2

    # here we filter on the minimum number of columns each taxa is in if requested
    if minrowsum > 0:
        df1_presence = df1_presence.loc[df1_presence[df1_presence.sum(axis=1) > minrowsum].index]
        df2_presence = df2_presence.loc[df2_presence[df2_presence.sum(axis=1) > (2 * minrowsum)].index]

    # combine the two matrices and sort them
    both = df1_presence.add(df2_presence, fill_value=0)
    sboth = both.loc[both.sum(axis=1).sort_values(ascending=asc_sort).index]
    sboth = sboth.sort_index(axis=1) # sort by column names

    return sboth

def read_the_data(sequence_type, datadir, sslevel='subsystems_norm_ss.tsv.gz', taxa="family", verbose=False):
    """
    Read the data and return the data frame and metadata
    :param sequence_type: MGI or MinION
    :param datadir: where is the data
    :param sslevel: subsystems level to read
    :param taxa: taxonomy level to read
    :return: two dataframes, data and metadata
    """

    ss_df = read_subsystems(
        os.path.join(datadir, sequence_type, "FunctionalAnalysis", "subsystems", sslevel), sequence_type)
    ss_df = ss_df.T
    if verbose:
        print(f"Read {ss_df.shape[0]} samples and {ss_df.shape[1]} subsystems", file=sys.stderr)
    genus_otu = read_taxonomy(datadir, sequence_type, taxa)
    genus_otu = genus_otu.T
    if verbose:
        print(f"Read {genus_otu.shape[0]} samples and {genus_otu.shape[1]} {taxa}", file=sys.stderr)
    df = ss_df.merge(genus_otu, left_index=True, right_index=True, how='inner')

    metadata = read_metadata(datadir, sequence_type, categorise=True)
    if verbose:
        print(f"Read {metadata.shape[0]} samples and {metadata.shape[1]} metadata columns", file=sys.stderr)

    return df, metadata
