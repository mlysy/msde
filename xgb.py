import xgboost as xgb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA, FastICA
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score

trainData = pd.read_csv("train1.csv")
test = pd.read_csv("test1.csv")

#trainSet, testSet = train_test_split(trainData,test_size=0.3)
#testY = testSet["y"]
#trainY = trainSet["y"]
trainY = trainData['y']

# Dimensionality Reduction******************************************
n_comp = 10

pca = PCA(n_components = n_comp)
#pca2_results_train = pca.fit_transform(trainSet.drop(["ID","y"], axis=1))
pca2_results_train = pca.fit_transform(trainData.drop(["ID","y"], axis=1))
#pca2_results_test = pca.transform(testSet.drop(["ID","y"], axis=1))
pca2_results_test = pca.transform(test.drop(["ID"], axis=1))

ica = FastICA(n_components = n_comp)
#ica2_results_train = ica.fit_transform(trainSet.drop(["ID","y"],axis=1))
ica2_results_train = ica.fit_transform(trainData.drop(["ID","y"],axis=1))
#ica2_results_test = ica.transform(testSet.drop(["ID","y"],axis=1))
ica2_results_test = ica.transform(test.drop(["ID"], axis=1))

for ii in range(1,n_comp+1):
#    trainSet['pca_' + str(ii)] = pca2_results_train[:,ii-1]
    trainData['pca_' + str(ii)] = pca2_results_train[:,ii-1]
#    testSet['pca_' + str(ii)] = pca2_results_test[:,ii-1]
#    trainSet['ica_' + str(ii)] = ica2_results_train[:,ii-1]
#    testSet['ica_' + str(ii)] = ica2_results_test[:,ii-1]
    test['pca_' + str(ii)] = pca2_results_test[:,ii-1]
    trainData['ica_' + str(ii)] = ica2_results_train[:,ii-1]
    test['ica_' + str(ii)] = ica2_results_test[:,ii-1]

#*******************************************************************
# Feature selection using Lasso*************************************

clf = Lasso(alpha = 0.1)

trainX = trainData.drop(['ID','y'],axis=1)
#trainX = trainSet.drop(['ID','y'],axis=1)
trainY = np.ravel(trainY)
#testY = np.ravel(testY)

clf.fit(trainX,trainY)

sigCoef = trainX.columns[abs(clf.coef_)>0.0]
#trainSet = trainSet[sigCoef]
trainSet = trainData[sigCoef]
testSet = test[sigCoef]
y_mean = np.mean(trainY)

xgb_params = {
        'n_trees': 500,
        'eta': 0.005,
        'max_depth': 4,
        'subsample':0.95,
        'objective':'reg:linear',
        'eval_metric': 'rmse',
        'base_score': y_mean,
        'silent':1
}

dtrain = xgb.DMatrix(trainSet,trainY)
dtest = xgb.DMatrix(testSet)

cv_result = xgb.cv(xgb_params,
                    dtrain,
                    num_boost_round=750,
                    early_stopping_rounds=50,
                    verbose_eval=50,
                    show_stdv=False)


num_boost_rounds = len(cv_result)
model = xgb.train(dict(xgb_params,silent=0),dtrain,num_boost_round=num_boost_rounds)

#print(r2_score(testY,model.predict(dtest)))
predY = model.predict(dtest)
submit = pd.DataFrame(test['ID'])
submit['y'] = pd.Series(predY)

submit.to_csv('submission.csv',index=False)
