import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef, balanced_accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold

def create_balanced_folds(y, n_splits, shuffle=True, random_state=None):
    """
    Create balanced folds for cross-validation by ensuring each fold has a proportional representation of class labels.

    Args:
        y (array-like): Array of class labels.
        n_splits (int): Number of cross-validation folds.
        shuffle (bool, optional): Whether to shuffle the data before creating folds. Defaults to True.
        random_state (int, optional): Seed for random number generator. Defaults to None.

    Returns:
        list: List of (train indices, test indices) tuples for each fold.

    """
    unique_labels, label_counts = np.unique(y, return_counts=True)
    folds = []
    
    for label, count in zip(unique_labels, label_counts):
        kf = KFold(n_splits=n_splits, shuffle=shuffle, random_state=random_state)
        label_indices = np.where(y == label)[0]
        label_folds = list(kf.split(label_indices))
        
        if len(folds) == 0:
            folds = label_folds
        else:
            for i in range(n_splits):
                train_indices, test_indices = label_folds[i]
                folds[i][0].extend(label_indices[train_indices])
                folds[i][1].extend(label_indices[test_indices])
    
    return(folds)

def custom_learning_curve(clf, X, y, train_sizes, cv):
    """
    Compute learning curve for a classifier using custom cross-validation.

    Args:
        clf (estimator): Classifier or estimator object.
        X (array-like): Input features.
        y (array-like): Target labels.
        train_sizes (array-like): List of training set sizes to evaluate.
        cv (int or cross-validation generator): Cross-validation strategy.

    Returns:
        dict: Dictionary containing learning curve results.

    """
    train_size_list = list()
    train_ra_scores, train_mc_scores, train_ba_scores  = list(), list(), list()
    valid_ra_scores, valid_mc_scores, valid_ba_scores  = list(), list(), list()
    
    train_sizes = [size for size in train_sizes if size <= len(list(cv.split(X, y))[0][0])]
    for train_size in train_sizes:

        for train_index, test_index in cv.split(X, y):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            # Train subset that respects class imbalance
            X_train_subset, _, y_train_subset, _ = train_test_split(X_train, y_train, train_size=train_size, stratify=y_train)
            # Train on the subset and evaluate the performance
            clf.fit(X_train_subset, y_train_subset)

            # Metrics for evaluating the model; account for class imbalance
            y_train_subset = y_train_subset.replace({'OVC':1,'Normal':0})
            y_pred_subset = pd.Series(clf.predict(X_train_subset)).replace({'OVC':1,'Normal':0})
            train_ra_score = roc_auc_score(y_train_subset, y_pred_subset)
            train_mc_score = matthews_corrcoef(y_train_subset, y_pred_subset)
            train_ba_score = balanced_accuracy_score(y_train_subset, y_pred_subset)

            y_test = y_test.replace({'OVC':1,'Normal':0})
            y_pred = pd.Series(clf.predict(X_test)).replace({'OVC':1,'Normal':0})
            valid_ra_score = roc