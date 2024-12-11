"""
Create some clean and encode data functions
"""
import re


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
