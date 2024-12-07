"""
Read some information about the metadata, like what the columns are, and what the values are
"""

import os
import sys
import gzip


__author__ = 'Rob Edwards'

def metadata_types(filename='Column Definitions.txt', metadata_path=os.path.join(os.path.dirname(__file__), '..', 'Metadata')):
    """
    Read the metadata file and return a data frame
    """
    metadata = os.path.join(metadata_path, filename)
    if not os.path.exists(metadata):
        print(f"Error: {metadata} does not exist", sys.stderr)
        return None
    print(f"Reading metadata from {metadata}", sys.stderr)

    types = {}

    with open(metadata, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            p = l.strip().split("\t")
            types[p[0]] = p[4]

    return types