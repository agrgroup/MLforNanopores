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

# Extract the targets (nanopore probabilities) and the nanopore features from the training set
y = train['probability']
X = train.drop(columns='probability')

# Extract the targets (nanopore probabilities) and the nanopore features from the test set
ytest = test['probability']
Xtest = test.drop(columns=['probability'])

# Plot log-histogram of probabilities
def plot_loghist(x, bins):
    _, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(x, bins=logbins, log=True)
    plt.xscale('log')

plot_loghist(ytest, 100)
plt.xlabel('True probability (log scale)', fontsize=15)
plt.ylabel('Frequency (log scale)', fontsize=15)
plt.savefig('plots/test_hist.png', dpi=1200, bbox_inches="tight")
plt.show()

# Load the pretrained model and test it using the test set
model = joblib.load('saved_model/probabilityPredictor.pkl')
test_preds = model.predict(Xtest)
test_r2 = metrics.r2_score(ytest, test_preds)

# Get the training predictions 
X = X[y<=0.4].reset_index(drop=True)
y = y[y<=0.4].reset_index(drop=True)
train_preds = model.predict(X)

# Print the r2 score obtained after testing the model on the test set
print(metrics.r2_score(ytest, test_preds))

# Scatter plot
plt.scatter(y, train_preds, color='purple', alpha=0.5, label='train')
plt.scatter(ytest, test_preds, color='green', alpha=0.5, label='test')
x = (np.min(ytest), np.max(ytest))
plt.plot(x, x, 'r', ls='dashed', label='y=x', linewidth=3.0)
plt.xlabel('True probability', fontsize=15)
plt.ylabel('Predicted probability', fontsize=15)
plt.annotate('R2'.translate(superscript) + f': {str(test_r2)[:6]}', xy=(80, 250), xycoords='axes points', 
            ha='right', va='top', bbox=dict(boxstyle='round', fc='w', ec='black', alpha=0.1), fontsize=15)
fig = plt.gcf()
plt.legend(loc='lower right', fontsize=15)
plt.show()
fig.savefig('plots/LGBM_tw_test_sp.png', bbox_inches="tight")