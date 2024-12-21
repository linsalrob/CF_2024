"""
The gradient boosting notebook, but as a script so we can run it for extended periods!
"""

import os
import sys
import argparse

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

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler, LabelEncoder, OrdinalEncoder
from sklearn.inspection import permutation_importance
from sklearn.impute import SimpleImputer

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.multioutput import MultiOutputClassifier, MultiOutputRegressor, ClassifierChain
from sklearn.metrics import mean_squared_error

from scipy.stats import linregress

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor


# there is a FutureWarning in sklearn StandardScalar which is really annoying. This ignores it.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from socket import gethostname
hostname = gethostname()
if hostname.startswith('hpc-node'):
    IN_DEEPTHOUGHT = True
    sys.path.append('..')
else:
    IN_DEEPTHOUGHT = False

import cf_analysis_lib

__author__ = 'Rob Edwards'


def gb_classifier(X, y, n_estimators=200):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = GradientBoostingClassifier(
        max_features="sqrt",
        n_estimators=n_estimators,
        learning_rate=0.005,
        min_samples_leaf=10,
        max_depth=5
    )

    """
    n_estimators=200,       # number of trees
    learning_rate=0.1,      # shrinkage
    max_depth=3,            # tree depth
    max_features='sqrt',    # random subset of features at each split
    min_samples_leaf=5,     # minimum samples per leaf
    random_state=42         # for reproducibility
    """

    model.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = model.predict(X_test)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)

    # Feature importance
    feature_importances = pd.DataFrame(model.feature_importances_, index=X.columns, columns=['importance'])
    feature_importances_sorted = feature_importances.sort_values(by='importance', ascending=False)
    return mse, feature_importances_sorted


def gb_regressor(X, y, n_estimators=200):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = GradientBoostingRegressor(max_features="sqrt", n_estimators=n_estimators)
    model.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = model.predict(X_test)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)

    # Feature importance
    feature_importances = pd.DataFrame(model.feature_importances_, index=X.columns, columns=['importance'])
    feature_importances_sorted = feature_importances.sort_values(by='importance', ascending=False)
    return mse, feature_importances_sorted


def random_forest_regression(X, y, n_estimators=200):
    """
    Run a regressor for continuous data and return the mean squared error and the feature importances
    """

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Initialize and train a RandomForestRegressor model
    model = RandomForestRegressor(random_state=42, n_estimators=n_estimators)  # You can adjust hyperparameters
    model.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = model.predict(X_test)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)

    # Feature importance
    feature_importances = pd.DataFrame(model.feature_importances_, index=X.columns, columns=['importance'])
    feature_importances_sorted = feature_importances.sort_values(by='importance', ascending=False)
    return mse, feature_importances_sorted


def random_forest_classifier(X, y, n_estimators=200):
    """
    Run a classifier for categorical data and return the mean squared error and the feature importances
    """

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Initialize and train a RandomForestRegressor model
    model = RandomForestClassifier(random_state=42, n_estimators=n_estimators)  # You can adjust hyperparameters
    model.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = model.predict(X_test)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)

    # Feature importance
    feature_importances = pd.DataFrame(model.feature_importances_, index=X.columns, columns=['importance'])
    feature_importances_sorted = feature_importances.sort_values(by='importance', ascending=False)
    return mse, feature_importances_sorted


def plot_feature_importance(ax, feature_importances_sorted, title):
    # Create dotted lines and circles for each feature
    for feature in feature_importances_sorted.index[::-1]:
        importance = feature_importances_sorted.loc[feature, 'importance']
        ax.plot([importance], [feature], linestyle='dotted', marker='o', markersize=5, c='blue')
        ax.plot([0, importance], [feature, feature], linestyle='dotted', marker='None', markersize=5, c='lightblue')

    ax.set_xlabel("Importance")
    ax.set_ylabel("")
    ax.set_title(title)


def plot_feature_abundance(ax, feature_df, intcol, title):
    """
    Plot the top n important features.

    use something like this:
    top20 = list(feature_importances_sorted[:20].index)+[intcol]
    plot_feature_abundance(ax, merged_df[top20], intcol, f"Plot of normalised measures that are important to distinguish '{intcol}' usage")
    """

    # before we plot the data we scale the data to make the mean 0 and the variance 1.
    # you can compare the values before and after by looking at merged_df[top20].max() and  scaled_df.max()
    # scaler = StandardScaler()
    scaler = MinMaxScaler()
    scaled_df = pd.DataFrame(scaler.fit_transform(feature_df), columns=feature_df.columns)
    scaled_df[intcol] = feature_df[intcol].values

    melted_df = pd.melt(scaled_df, id_vars=[intcol], var_name='Feature', value_name='Value')

    sns.boxplot(data=melted_df, x='Value', y='Feature', hue=intcol, fill=False, legend=False, color='k', fliersize=0,
                ax=ax)
    sns.stripplot(data=melted_df, x='Value', y='Feature', hue=intcol, jitter=True, alpha=0.5, dodge=True, ax=ax)

    ax.set_title(title)
    ax.set_xlabel('Normalised Abundance')
    ax.set_ylabel('')

def read_the_data(sequence_type, datadir, sslevel = 'subsystems_norm_ss.tsv.gz', taxa = "family"):

    ss_df = cf_analysis_lib.read_subsystems(
        os.path.join(datadir, sequence_type, "FunctionalAnalysis", "subsystems", sslevel), sequence_type)
    ss_df = ss_df.T
    genus_otu = cf_analysis_lib.read_taxonomy(datadir, sequence_type, taxa)
    genus_otu = genus_otu.T
    df = ss_df.merge(genus_otu, left_index=True, right_index=True, how='inner')

    metadata = cf_analysis_lib.read_metadata(datadir, sequence_type, categorise=True)

    return df, metadata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--results', help='results file', required=True)
    parser.add_argument('-n', '--iterations', help='number of iterations', type=int, default=100)
    parser.add_argument('-i', '--images', help='image directory', default='gb_images')
    parser.add_argument('-d', '--datadir', help='data directory', default='..')
    parser.add_argument('-s', '--sequence_type', help='sequence type', default='MGI')
    parser.add_argument('-l', '--sslevel', help='subsystem level', default='subsystems_norm_ss.tsv.gz')
    parser.add_argument('-t', '--taxa', help='taxonomic level', default='family')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if os.path.exists(args.r):
        print(f"Error: {args.r} already exists. Cowardly not overwriting", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.images, exist_ok=True)

    df, metadata = read_the_data(sequence_type=args.sequence_type, datadir=args.datadir, sslevel=args.sslevel,
                                 taxa=args.taxa)

    replace_index = re.compile(r'^\d+\s+')
    replace_nonword = re.compile(r'\W+')

    resultsfile = open(args.r, 'w')
    print(f"Predictor\tFeature\tAvg. Importance\tNumber of iterations (out of {args.iterations})", file=resultsfile)

    for intcol in df.columns:
        print(f"Working on {intcol}", file=sys.stderr)
        intcol_title = replace_index.sub('', intcol)
        intcol_filename = intcol.replace(" ", "_")
        intcol_filename = replace_nonword.sub('', intcol_filename)

        merged_df = df.join(metadata[[intcol]])

        X = merged_df.drop(intcol, axis=1)
        y = merged_df[intcol]

        top_features = {}
        top_feature_counts = {}
        for i in range(args.iterations):
            if metadata[intcol].dtype == 'object':
                mse, feature_importances_sorted = gb_classifier(X, y)
                met = 'classifier'
            else:
                mse, feature_importances_sorted = gb_regressor(X, y)
                met = 'regressor'

            for f in feature_importances_sorted.index[:n]:
                top_features[f] = top_features.get(f, 0) + feature_importances_sorted.loc[f, 'importance']
                top_feature_counts[f] = top_feature_counts.get(f, 0) + 1

        sorted_top_feats = sorted(top_features, key=lambda x: top_features[x], reverse=True)
        for x in sorted_top_feats[:10]:
            print(f"{intcol}\t{x}\t{top_features[x] / args.iterations : .4f}\t\t{top_feature_counts[x]}", file=resultsfile)


        tfdf = pd.DataFrame.from_dict(top_features, orient="index", columns=["importance"]).sort_values(by='importance',
                                                                                                        ascending=False)
        n = 20
        topN = list(tfdf[:n].index) + [intcol]
        fig, axes = plt.subplots(figsize=(10, 6), nrows=1, ncols=2, sharey='row', sharex='col')
        plot_feature_importance(axes[0], tfdf[:n][::-1], "")
        plot_feature_abundance(axes[1], merged_df[topN][::-1], intcol, intcol_title)

        custom_labels = {0: 'No', 1: 'Yes'}
        handles, labels = axes[1].get_legend_handles_labels()  # Get one set of handles and labels
        updated_labels = [custom_labels[int(label)] for label in labels]

        for ax in axes.flat:
            if ax.get_legend() is not None:  # Check if legend exists
                ax.get_legend().remove()

        plt.xticks(rotation=90)
        fig.legend(handles, updated_labels, loc='upper center', ncol=2, fontsize=12, title="Antibiotic Usage")
        plt.tight_layout(rect=[0, 0, 1, 0.9])
        plt.savefig(f"{args.images}/{intcol_filename}.png")

    resultsfile.close()