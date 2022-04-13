import pandas as pd
import numpy as np

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

# Load the training and test sets
data = pd.read_csv('input/train.csv')
test = pd.read_csv('input/test.csv')

# Extract the targets (nanopore probability) and features (nanopores descriptors) from the training set
y = data['probability']
X = data.drop(columns=['probability'])

# Plot log-histogram of y
def plot_loghist(x, bins):
    _, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(x, bins=logbins, log=True)
    plt.xscale('log')

plot_loghist(y, 100)
plt.xlabel('True probability (log scale)', fontsize=15)
plt.ylabel('Frequency (log scale)', fontsize=15)
plt.savefig('plots/hist.png', dpi=1200, bbox_inches="tight")
plt.show()

# Tune hyperparamaters
# model = Model()
# catboost_params, lgbm_params = model.tune(X, y, iter=10)

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

# Train the two-stage model and save it
model = Model(lgbm_params, catboost_params)
model.fit(X, y)
model.save(filename='saved_model/probabilityPredictor.pkl')
