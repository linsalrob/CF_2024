"""
Read some information about the metadata, like what the columns are, and what the values are
"""

import os
import sys


__author__ = 'Rob Edwards'

def metadata_definitions(filename='Column Definitions.txt', metadata_path=os.path.join(os.path.dirname(__file__), '..', 'Metadata')):
    """
    Read the metadata file and return a data frame
    """
    metadata_file = os.path.join(metadata_path, filename)
    if not os.path.exists(metadata_file):
        print(f"Error: {metadata_file} does not exist", sys.stderr)
        return None

    types = {}

    with open(metadata_file, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            p = l.strip().split("\t")
            types[p[1]] = p[4]

    return types
