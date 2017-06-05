import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder

def clean_data(raw_data):
    # Initialising empty data frame for output
    clean = pd.DataFrame()

    # Loop through columns to replace categorical data with labels
    for column in raw_data:
        current_column = raw_data[column]
        if current_column.dtype == "object":
            lbl = LabelEncoder()
            lbl.fit(list(current_column)) # Fit the labels of current column
            clean[column] = lbl.transform(list(current_column)) # return labels in from 0,...,n-1
    
    return clean

#trainData = pd.read_csv("train.csv")
#testData = pd.read_csv("test.csv")
#
#ipdb.set_trace()
#
#for colIndex in trainData.columns:
#    if trainData[colIndex].dtype == 'object':
#        lbl.fit(list(trainData[colIndex].values) + list(testData[colIndex].values))
#        trainData[colIndex] = lbl.transform(list(trainData[colIndex].values))
#        testData[colIndex] = lbl.transform(list(testData[colIndex].values))
#
#trainData.to_csv('train1.csv', index = False)
#testData.to_csv('test1.csv', index = False)

#for colIndex in xrange(2,10):
#    uniqueId = trainData.iloc[:,colIndex].unique()
#    tmpLen = len(uniqueId)
#    for rowIndex in xrange(0,tmpLen):
#       addC = trainData.iloc[:,colIndex] == uniqueId[rowIndex]
#       addC = addC.astype(int)
#       aname = 'X' + str(colIndex) + '-' + str(rowIndex)
#       trainData[aname] = pd.Series(addC,index=trainData.index)
#

#newTrain = trainData.iloc[:,10:]
#newTrain["ID"] = pd.Series(trainData.iloc[:,0],index=trainData.index)
#newTrain["y"] = pd.Series(trainData.iloc[:,1],index=trainData.index)


#newTrain.to_csv('train1.csv', index = False)
