"""
Random Forests and Gradient Boosted Random Forests for classification and regression

This was cut from the GBRF jupyter notebook so we can use it elsewhere.

"""

import pandas as pd

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

