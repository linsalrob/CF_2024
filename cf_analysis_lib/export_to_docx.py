"""
Export a pandas data frame to MS Word so we can use
it in the paper!
"""

import sys
from docx import Document

def pd2docx(df, docx, verbose=False):
    """
    Export a pandas data frame to a docx file

    :param df: the pandas data frame
    :param docx: the name of the docx file
    :param verbose: more output
    :return: None
    """

    if verbose:
        print(f"Writing the data frame to {docx}", file=sys.stderr)

    # Create a new Word Document
    document = Document()

    # Add a table with one extra row for the header
    table = document.add_table(rows=1, cols=len(df.columns))
    header_cells = table.rows[0].cells
    for i, column in enumerate(df.columns):
        header_cells[i].text = str(column)

    # Add the DataFrame rows to the table
    for index, row in df.iterrows():
        row_cells = table.add_row().cells
        for i, cell in enumerate(row):
            row_cells[i].text = str(cell)

    # Save the document
    document.save(docx)
    return None
