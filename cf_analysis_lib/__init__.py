"""
A library of functions for analyzing the output of the CF analysis pipeline.
"""


import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
import json
import random

from itertools import cycle
from collections import Counter

import networkx as nx
import community

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler, LabelEncoder, OrdinalEncoder
from sklearn.inspection import permutation_importance
from sklearn.impute import SimpleImputer

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.multioutput import MultiOutputClassifier, MultiOutputRegressor, ClassifierChain
from sklearn.metrics import mean_squared_error, roc_curve, auc

from scipy.stats import linregress
from scipy.cluster.hierarchy import linkage, fcluster

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor


# there is a FutureWarning in sklearn StandardScalar which is really annoying. This ignores it.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

try:
  import google.colab
  IN_COLAB = True
  from google.colab import drive
  drive.mount('/content/drive')
  datadir = '/content/drive/MyDrive/Projects/CF/Adelaide/CF_Data_Analysis'
except ImportError:
  IN_COLAB = False
  datadir = '..'


from socket import gethostname
hostname = gethostname()
if hostname.startswith('hpc-node'):
    IN_DEEPTHOUGHT = True
    sys.path.append('..')
else:
    IN_DEEPTHOUGHT = False


from adjustText import adjust_text




from .read_data import read_taxonomy, read_metadata, corrections, read_subsystems, sorted_presence_absence, read_the_data
from .metadata_data import metadata_definitions
from .clean_data import compatible_columns, categories_to_numeric, remove_highly_correlated_data
from .pathogens import BacterialPathogens
from .random_forests import random_forest_classifier, random_forest_regression, gb_classifier, gb_regressor
from .plot_figures import plot_pca, plot_feature_importance, plot_feature_abundance, plot_top_features
from .plot_figures import plot_abundance_stripplot, create_custom_labels
from .export_to_docx import pd2docx

__all__ = ['read_taxonomy', 'read_metadata', 'corrections', 'read_subsystems',
           'sorted_presence_absence', 'compatible_columns', 'categories_to_numeric', 'remove_highly_correlated_data',
           'metadata_definitions', 'read_the_data',
           'BacterialPathogens',
           'random_forest_classifier', 'random_forest_regression', 'gb_classifier', 'gb_regressor', 'plot_top_features',
           'plot_pca', 'plot_feature_importance', 'plot_feature_abundance', 'plot_top_features', 'plot_abundance_stripplot',
           'pd2docx', 'create_custom_labels'
        ]
