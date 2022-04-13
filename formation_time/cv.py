import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn import metrics

from model import Model

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

plt.rc('font', family='Arial')
plt.rc('font', size=25)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)

import warnings
warnings.filterwarnings('ignore')

# Threshold to split the dataset for training in two stages
t = 140

# Load the training set
train = pd.read_csv('input/train.csv')

# Extract the targets (nanopore formation times) and features (nanopores descriptors)
X, y = train.drop(columns='time'), train['time']

# Tuned parameters
catboost_params = {
    'bagging_temperature': 0.5702320183528147,
    'border_count': 62,
    'depth': 1,
    'l2_leaf_reg': 20,
    'learning_rate': 0.997203113617485,
    'od_type': 'IncToDec',
    'od_wait': 34,
    'random_strength': 74.14613739817753,
    'scale_pos_weight': 85.03369763670949
}

lgbm_params = {
    'colsample_bytree': 0.8042072656392166,
    'learning_rate': 0.6442221474121814,
    'max_bin': 758,
    'max_depth': 24,
    'num_leaves': 56,
    'reg_alpha': 0.04508502532105237,
    'reg_lambda': 19.36777109067546,
    'subsample': 0.6837522522496777,
}

classes = [1 if y_>t else 0 for y_ in y]

# Cross validation
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=40)
r2s, maes = [], []

for fold, (train_idx, valid_idx) in enumerate(skf.split(X, classes)):
    print(f'Fold: {fold}')

    # Split the dataset into training and validation sets
    Xtrain, Xvalid = X.iloc[train_idx].reset_index(drop=True), X.iloc[valid_idx].reset_index(drop=True)
    ytrain, yvalid = y.iloc[train_idx].reset_index(drop=True), y.iloc[valid_idx].reset_index(drop=True)
    
    # Initialize the model and train it
    model = Model(lgbm_params=lgbm_params, catboost_params=catboost_params)
    model.fit(Xtrain, ytrain)
    
    # Use the model to predict on the test set
    yval_pred = model.predict(Xvalid)
    
    # Calculate metrics such as r2 and mae 
    r2 = metrics.r2_score(yvalid, yval_pred)
    print(r2)
    r2s.append(r2)

    mae = metrics.mean_absolute_error(yvalid, yval_pred)
    print(mae)
    maes.append(mae)

# Calculate the mean r2 and mae over all folds
print(f'Overall_R2_score: {np.mean(r2s)}, Overall_MAE_score: {np.mean(maes)}')