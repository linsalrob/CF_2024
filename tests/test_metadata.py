# tests/test_metadata.py
import os
import unittest
from cf_analysis_lib.metadata_data import metadata_definitions

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.metadata_path = os.path.join(os.path.dirname(__file__), '..', 'Metadata')
        self.filename = 'Column Definitions.txt'

    def test_metadata_types(self):
        types = metadata_definitions(self.filename, self.metadata_path)
        self.assertIsNotNone(types)
        self.assertIsInstance(types, dict)

if __name__ == '__main__':
    unittest.main()