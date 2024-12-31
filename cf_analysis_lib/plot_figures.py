"""
Plot a PCA and a read abundance plot for the interesting clusters
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns

from sklearn.preprocessing import MinMaxScaler, OrdinalEncoder

from itertools import cycle
from sklearn.decomposition import PCA
from adjustText import adjust_text

__author__ = 'Rob Edwards'


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
    tmpdf = feature_df.drop(intcol, axis=1)
    scaled_df = pd.DataFrame(scaler.fit_transform(tmpdf), columns=tmpdf.columns)
    scaled_df[intcol] = feature_df[intcol].values

    melted_df = pd.melt(scaled_df, id_vars=[intcol], var_name='Feature', value_name='Value')

    sns.boxplot(data=melted_df, x='Value', y='Feature', hue=intcol, fill=False, legend=False, color='k', fliersize=0,
                ax=ax)
    sns.stripplot(data=melted_df, x='Value', y='Feature', hue=intcol, jitter=True, alpha=0.5, dodge=True, ax=ax)

    ax.set_title(title)
    ax.set_xlabel('Normalised Abundance')
    ax.set_ylabel('')


def plot_top_features(merged_df, top_features, top_feature_counts, intcol, intcol_title, custom_labels=None):
    sorted_top_feats = sorted(top_features, key=lambda x: top_features[x], reverse=True)
    for x in sorted_top_feats[:10]:
        print(f"{x}: {top_features[x] / 5:.4f} ({top_feature_counts[x]} times)")

    n = 20
    tfdf = pd.DataFrame.from_dict(top_features, orient="index", columns=["importance"]).sort_values(by='importance',
                                                                                                    ascending=False)
    topN = list(tfdf[:n].index) + [intcol]
    fig, axes = plt.subplots(figsize=(10, 6), nrows=1, ncols=2, sharey='row', sharex='col')
    plot_feature_importance(axes[0], tfdf[:n][::-1], "")
    plot_feature_abundance(axes[1], merged_df[topN][::-1], intcol, intcol_title)

    if not custom_labels:
        custom_labels = {0: 'No', 1: 'Yes'}

    handles, labels = axes[1].get_legend_handles_labels()  # Get one set of handles and labels
    updated_labels = [custom_labels[float(label)] for label in labels]

    for ax in axes.flat:
        if ax.get_legend() is not None:  # Check if legend exists
            ax.get_legend().remove()

    plt.xticks(rotation=90)
    fig.legend(handles, updated_labels, loc='upper center', ncol=2, title=intcol_title)
    plt.tight_layout(rect=[0, 0, 1, 0.9])
    plt.show()

def plot_pca(ax, df, metadata, cluster_assignments, interesting_cluster, intcol):
    pca = PCA(n_components=2)

    df_clust = df[cluster_assignments.loc[cluster_assignments["Cluster"] == interesting_cluster, "Feature"]]
    merged_df_clust = df_clust.join(metadata[[intcol]])

    pca_result = pca.fit_transform(df_clust)
    pca_df = pd.DataFrame(data=pca_result, index=df.index, columns=['PC1', 'PC2'])

    # Get loadings
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loadings_df = pd.DataFrame(loadings, index=df_clust.columns, columns=['PC1', 'PC2'])

    # Create a DataFrame for top loadings
    top_loadings_df = loadings_df.loc[loadings_df['PC1'].abs().sort_values(ascending=False).index]
    top_loadings_df.head()

    explained_variance = pca.explained_variance_ratio_ * 100
    pc1_variance = explained_variance[0]
    pc2_variance = explained_variance[1]

    sns.scatterplot(data=pca_df, x='PC1', y='PC2', alpha=0.2, ax=ax)
    ax.set_title(f"Cluster {interesting_cluster}")
    ax.set_xlabel(f'Principal Component 1 ({pc1_variance:.3f}%)')
    ax.set_ylabel(f'Principal Component 2 ({pc2_variance:.3f}%)')

    # add the loadings ... we only plot maxloadings here
    maxloadings = 5
    if len(loadings) < maxloadings:
        maxloadings = len(loadings)

    plotscaler = 2
    texts = []
    colour_cycle = cycle(mcolors.TABLEAU_COLORS)
    found_pseudomonas = False
    for i in range(maxloadings):
        c = next(colour_cycle)
        xpos = top_loadings_df.iloc[i, 0] * plotscaler
        ypos = top_loadings_df.iloc[i, 1] * plotscaler
        ax.arrow(0, 0, xpos, ypos,
                      color=c, alpha=0.5, width=0.05)
        texts.append(ax.text(xpos, ypos, top_loadings_df.index[i], color=c))

    adjust_text(texts, ax=ax)

def plot_abundance_stripplot(ax, df, metadata, cluster_assignments, interesting_cluster, intcol):
    df_clust = df[cluster_assignments.loc[cluster_assignments["Cluster"] == interesting_cluster, "Feature"]]
    merged_df_clust = df_clust.join(metadata[[intcol]])
    df_clust_m = merged_df_clust.melt(id_vars=intcol, var_name='Features', value_name='Normalised read abundance')
    sns.stripplot(data=df_clust_m, x='Features', y='Normalised read abundance',
                  hue=intcol, dodge=True, jitter=True, ax=ax)
    ax.tick_params(axis='x', rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment('right')
    ax.set_title(f"Read abundance for {intcol} cluster {interesting_cluster}")
