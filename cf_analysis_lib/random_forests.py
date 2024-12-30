"""
Random Forests and Gradient Boosted Random Forests for classification and regression

This was cut from the GBRF jupyter notebook so we can use it elsewhere.

"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.preprocessing import MinMaxScaler, OrdinalEncoder

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier, GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error





def gb_classifier(X, y, n_estimators=10000):
    """
    Run a classifier for categorical data and return the mean squared error and the feature importances
    """

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


def gb_regressor(X, y, n_estimators=10000):
    """
    Run a regressor for continuous data and return the mean squared error and the feature importances
    """

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


def random_forest_regression(X, y, n_estimators=10000):
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


def random_forest_classifier(X, y, n_estimators=10000):
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
