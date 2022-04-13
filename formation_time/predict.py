import pandas as pd
import numpy as np
import joblib

from sklearn import metrics

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

plt.rc('font', family='Arial')
plt.rc('font', size=25)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)

import warnings
warnings.filterwarnings('ignore')

superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

# Load the datasets
train = pd.read_csv('input/train.csv')
test = pd.read_csv('input/test.csv')

# Extract the targets (nanopore formation times) and the nanopore features from the training set
y = train['time']
X = train.drop(columns='time')

# Extract the targets (nanopore formation times) and the nanopore features from the test set
ytest = test['time']
Xtest = test.drop(columns=['time'])

# Plot histogram of probabilities
plt.hist(test['time'], bins=100)
plt.xlabel('True formation time (s)', fontsize=15)
plt.ylabel('Frequency', fontsize=15)
plt.savefig('plots/test_hist.png', dpi=1200, bbox_inches="tight")
plt.show()

# Load the pretrained model and test it using the test set
model = joblib.load('saved_model/fomationTimePredictor.pkl')
test_preds = model.predict(Xtest)
test_r2 = metrics.r2_score(ytest, test_preds)

# Get the train predictions 
train_preds = model.predict(X)

# Print the r2 score obtained after testing the model on the test set
print(metrics.mean_absolute_error(y, train_preds), metrics.mean_absolute_error(ytest, test_preds))

# Scatter plot
plt.scatter(y, train_preds, color='purple', alpha=0.5, label='train')
plt.scatter(ytest, test_preds, color='green', alpha=0.5, label='test')
plt.xticks(np.arange(75, 250, step=25))
plt.yticks(np.arange(75, 250, step=25))
x = (np.min(ytest), np.max(ytest))
plt.plot(x, x, 'r', ls='dashed', label='y=x', linewidth=3.0)
plt.xlabel('True formation time (s)', fontsize=15)
plt.ylabel('Predicted formation time (s)', fontsize=15)
plt.annotate('R2'.translate(superscript) + f': {str(test_r2)[:6]}', xy=(80, 250), xycoords='axes points', 
            ha='right', va='top', bbox=dict(boxstyle='round', fc='w', ec='black', alpha=0.1), fontsize=15)
fig = plt.gcf()
plt.legend(loc='lower right', fontsize=15)
plt.show()
fig.savefig('plots/LGBM_rmse_test_sp.png', bbox_inches="tight")