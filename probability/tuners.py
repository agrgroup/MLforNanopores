import numpy as np

from sklearn.model_selection import KFold, StratifiedKFold
from sklearn import metrics

import optuna
from functools import partial

import lightgbm as lgb
import catboost as cb

def to_class(x, t=0.5):
    '''
    Function to convert the predicted probabilities into classes (0, 1) based on a threshold t
    '''
    x = 1/(1+np.exp(-x))
    y = np.zeros(len(x))
    y[x>=t] = 1
    y = y.astype('int')
    return y

def optimize_cb_classifier(trial, X, y):
    '''
    Function to tune the CatBoost classification model
    '''
    depth = trial.suggest_int('depth', 1, 8)
    learning_rate = trial.suggest_uniform('learning_rate', 1e-4, 1)
    random_strength = trial.suggest_uniform('random_strength', 1e-9, 100)
    bagging_temperature = trial.suggest_uniform('bagging_temperature', 1e-2, 0.99)
    border_count = trial.suggest_int('border_count', 1, 255)
    l2_leaf_reg = trial.suggest_int('l2_leaf_reg', 2, 200)
    scale_pos_weight = trial.suggest_uniform('scale_pos_weight', 1e-3, 100)
    od_type = trial.suggest_categorical('od_type', ['IncToDec', 'Iter'])
    od_wait = trial.suggest_int('od_wait', 10, 300)

    
    params = {
        'random_state': 42,
        'objective': 'Logloss',
        'random_strength': random_strength,
        'border_count': border_count,
        'learning_rate': learning_rate,
        'depth': depth,
        'scale_pos_weight': scale_pos_weight,
        'od_type': od_type,
        'od_wait': od_wait,
        'l2_leaf_reg': l2_leaf_reg,
        'bagging_temperature': bagging_temperature
    }
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    acc = []
    for train_idx, valid_idx in skf.split(X, y):
        
        X_train, X_valid = X[train_idx], X[valid_idx]
        y_train, y_valid = y[train_idx], y[valid_idx]
        
        PoolTrain = cb.Pool(data=X_train, label=y_train)
        PoolValid = cb.Pool(data=X_valid, label=y_valid)
    
        catboost = cb.train(
            params=params,
            pool=PoolTrain,
            num_boost_round=15000,  
            early_stopping_rounds=50,
            evals=(PoolValid),
            verbose=False
        )
        
        preds = to_class(catboost.predict(PoolValid))
        acc.append(metrics.balanced_accuracy_score(y_valid, preds))
    
    return np.mean(acc)

def optimize_lgbm_tweedie(trial, X, y):  
    '''
    Function to tune the LightGBM regression model
    '''  
    max_depth = trial.suggest_int('max_depth', -1, 1200)
    learning_rate = trial.suggest_uniform('learning_rate', 1e-5, 0.99)
    num_leaves = trial.suggest_int('num_leaves', 2, 300)
    subsample = trial.suggest_uniform('subsample', 0.01, 1)
    colsample_bytree = trial.suggest_uniform('colsample_bytree', 0.01, 1)
    reg_lambda = trial.suggest_uniform('reg_lambda', 0.01, 100)
    reg_alpha = trial.suggest_uniform('reg_alpha', 0.01, 100)
    max_bin = trial.suggest_int('max_bin', 50, 800)
    
    
    params = {
        'metric': 'tweedie',
        'objective': 'tweedie',
        'boosting_type': 'gbdt',
        'nthreads': -1,
        'num_leaves': num_leaves,
        'max_bin': max_bin, 
        'reg_alpha': reg_alpha, 
        'colsample_bytree': colsample_bytree, 
        'reg_lambda': reg_lambda, 
        'subsample': subsample, 
        'max_depth': max_depth, 
        'learning_rate': learning_rate, 
        'verbosity': 0, 
        'random_state': 42, 
    }
    
    oof_preds = np.zeros(len(X))
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    
    for train_idx, valid_idx in kf.split(X, y):
        
        X_train, X_valid = X.iloc[train_idx], X.iloc[valid_idx]
        y_train, y_valid = y.iloc[train_idx], y.iloc[valid_idx]
                
        LGBMTrain = lgb.Dataset(data=X_train, label=y_train)
        LGBMValid = lgb.Dataset(data=X_valid, label=y_valid)
    
        lgbm = lgb.train(
            params,
            LGBMTrain,
            num_boost_round=15000,
            valid_sets=[LGBMValid],
            early_stopping_rounds=50, 
            verbose_eval=False    
        )
        
        preds = lgbm.predict(X_valid) 
        oof_preds[valid_idx] = preds
    
    r2 = metrics.r2_score(y, oof_preds)    
    return r2

def tune_cb(X, y_c, iter=100):
    ''' 
    Wrapper function to call the "optimize_cb_classifier" function for tuning the CatBoost classification model
    '''
    optimization_func = partial(optimize_cb_classifier, X=X.values, y=y_c)
    study = optuna.create_study(direction='maximize', study_name='best_params', storage='sqlite:///optim_xgb.db', load_if_exists=True)
    study.optimize(optimization_func, n_trials=iter, n_jobs=1)
    return study.best_params

def tune_lgbm(X, y, t, iter=100):
    '''
    Wrapper function to call the "optimize_lgbm_tweedie" function for tuning the LightGBM regression model
    '''
    optimization_func = partial(optimize_lgbm_tweedie, X=X[y<t], y=y[y<t])
    study = optuna.create_study(direction='maximize', study_name='best_params', storage='sqlite:///optim_lgbm.db', load_if_exists=True)
    study.optimize(optimization_func, n_trials=iter, n_jobs=-1)
    return study.best_params

