"""
Create some clean and encode data functions
"""
import re
import pandas as pd


def compatible_columns(df):
    """
    Make really simple column names that do not have spaces and do not start with numbers
    :param df: The data frame
    :return:
    """

    # remove spaces from the column names
    df.columns = [re.sub(r'\s+', '_', c) for c in df.columns]
    # remove numbers from the start of the column names
    nos = {"1": "one", "2": "two", "3": "three", "4": "four", "5": "five", "6": "six", "7": "seven", "8": "eight",
           "9": "nine"}
    df.columns = [re.sub(r'^(\d)', lambda x: nos[x.group(1)] + "_", c) for c in df.columns]
    return df

def categories_to_numeric(df, sequence_type):
    """
    Convert categories and objects to numeric forms
    We use numbers if we have them, we keep the MGI/MinION ID so we can use it to join tables, but we don't want to use it as a predictor.
    We convert categories to their codes so we get numbers.
    We drop the other columns so we don't try and predict on free text columns.
    I'm not sure what to do about dates, yet, so I ignore them for now
    :param df: The data frame to clean
    :param sequence_type: The name of the column that is the sequence type
    :return: A clean data frame
    """

    # we make a copy of each data frame and then overwrite the columns. Internally this is more efficient than making a df one
    # column at a time.
    encoded = df.copy()
    to_delete = ['Pseudomonas_Culture']  # this is a duplicate columns with 'CS_Pseudomonas_aeruginosa'
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            continue
        elif col == sequence_type:
            continue
        elif df[col].dtypes == 'category':
            encoded[col] = df[col].cat.codes
        else:
            print(f'Dropping {col}: {df[col].dtypes}')
            to_delete.append(col)

    encoded = encoded.drop(columns=to_delete)
    return encoded
