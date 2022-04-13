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

# Extract the targets (nanopore formation times) and features (nanopores descriptors) from the training set
y = data['time']
X = data.drop(columns=['time'])

# Plot histogram of targets (nanopore formation times)
plt.hist(y, bins=100)
plt.xlabel('True formation time (s)', fontsize=15)
plt.ylabel('Frequency', fontsize=15)
plt.savefig('plots/hist.png', dpi=1200, bbox_inches="tight")
plt.show()

# # Tune hyperparamaters
# # model = Predictor()
# # catboost_params, lgbm_params = model.tune(X, y, iter=10)

# Tuned parameters
catboost_params = {
    'bagging_temperature': 0.3689921407560029,
    'border_count': 142,
    'depth': 7,
    'l2_leaf_reg': 40,
    'learning_rate': 0.45184710488781676,
    'od_type': 'IncToDec',
    'od_wait': 141,
    'random_strength': 54.85159387524285,
    'scale_pos_weight': 46.37952641283065
}

lgbm_params = {
    'colsample_bytree': 0.8347796498785935,
    'learning_rate': 0.06883504812393584,
    'max_bin': 95,
    'max_depth': 532,
    'num_leaves': 4,
    'reg_alpha': 43.9032829810609,
    'reg_lambda': 94.44388158539495,
    'subsample': 0.4568016573760509
}

# Train the two-stage model and save it
model = Model(lgbm_params, catboost_params, t=140)
model.fit(X, y)
model.save(filename='saved_model/fomationTimePredictor.pkl')
