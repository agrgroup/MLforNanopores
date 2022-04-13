import numpy as np
import joblib

import tuners

from sklearn.linear_model import LinearRegression

import lightgbm as lgb
import catboost as cb

class Model:
    '''
    2-stage model proposed to predict graphene nanopore probabilities
    '''
    def __init__(self, lgbm_params=dict(), catboost_params=dict(), t=0.015):
        '''
        lgbm_params: tuned lgbm (tweedie) params, which can otherwise be obtained using the "tune" attribute
        catboost_params: tuned catboost (classifier) params, which can otherwise be obtained using the "tune" attribute
        t: probability threshold for dealing (training and testing) seperately. 0.015 emperically worked the best
        '''
        self.t = t 
        self.lr = LinearRegression()
        self.lgbm = lgb.LGBMRegressor(
            random_state=42,
            n_estimators=15000,
            objective='tweedie',
            **lgbm_params
        )
        self.cb = cb.CatBoostClassifier(
            random_state=42,
            n_estimators=15000,
            verbose=False,
            **catboost_params,
        )
        
    def tune(self, X, y, iter=100):
        '''
        X: input features as a dataframe
        y: input probabilities
        iter: number of iteration until the tuning algorithm is required to run
        '''
        y_c = np.zeros(len(y))
        y_c[y >= self.t] = 1
        
        print('Tuning CatBoost....')
        tuned_cb_params = tuners.tune_cb(X, y_c, iter)
        
        print('Tuning LightGBM....')
        tuned_lgbm_params = tuners.tune_lgbm(X, y, self.t, iter)
        
        return tuned_cb_params, tuned_lgbm_params

    def stage_1(self, X, y):
        '''
        X: input features as a dataframe
        y: input probabilities  
        '''
        y_c = np.zeros(len(y))
        y_c[y>=self.t] = 1
        
        self.cb.fit(X, y_c)  

    def stage_2(self, X, y):
        '''
        X: input features as a dataframe
        y: input probabilities  
        '''
        self.lr.fit(X[y>=self.t], np.log(y[y>=self.t]))
        self.lgbm.fit(X[y<self.t], y[y<self.t])
        
    def fit(self, X, y):   
        '''
        X: input features as a dataframe
        y: input probabilities  
        '''         
        
        self.stage_1(X, y)
        self.stage_2(X, y)
                
    def save(self, filename):
        '''
        filename: filename for saving the model 
        '''
        joblib.dump(self, filename)
        
    def predict(self, Xtest):
        '''
        Xtest: test set features for which the probabilities are required
        '''
        
        binary_t = self.cb.predict(Xtest)
        
        Xtest_greater_t = Xtest[binary_t==1]
        Xtest_lesser_t = Xtest[binary_t==0]       
        
        ytest = np.zeros(len(Xtest))
        
        ytest[Xtest_greater_t.index] = np.exp(self.lr.predict(Xtest_greater_t))
        ytest[Xtest_lesser_t.index] = self.lgbm.predict(Xtest_lesser_t)
        
        return ytest