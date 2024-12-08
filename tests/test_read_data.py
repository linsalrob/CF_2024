# tests/test_read_data.py
import os
import sys
import unittest
import pandas as pd
from cf_analysis_lib.read_data import read_taxonomy, read_metadata, sorted_presence_absence, read_subsystems

class TestReadData(unittest.TestCase):

    def setUp(self):
        self.datadir = '.'
        self.sequence_type = 'MGI'
        self.taxonomy = 'genus'
        self.subsystems_file = os.path.join(self.datadir, self.sequence_type, 'FunctionalAnalysis', 'subsystems', 'subsystems_norm_ss.tsv.gz')

    def test_read_taxonomy(self):
        df = read_taxonomy(self.datadir, self.sequence_type, self.taxonomy)
        self.assertIsNotNone(df)
        self.assertIsInstance(df, pd.DataFrame)

    def test_read_metadata(self):
        df = read_metadata(self.datadir, self.sequence_type)
        self.assertIsNotNone(df)
        self.assertIsInstance(df, pd.DataFrame)

    def test_sorted_presence_absence(self):
        df1 = pd.DataFrame({'A': [1, 2], 'B': [3, 4]}, index=['tax1', 'tax2'])
        df2 = pd.DataFrame({'A': [0, 1], 'B': [1, 0]}, index=['tax1', 'tax2'])
        result = sorted_presence_absence(df1, df2)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)

    def test_read_subsystems(self):
        df = read_subsystems(self.subsystems_file, self.sequence_type)
        self.assertIsNotNone(df)
        self.assertIsInstance(df, pd.DataFrame)

if __name__ == '__main__':
    unittest.main()