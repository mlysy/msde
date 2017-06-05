import numpy as np
import pandas as pd
import sklearn as sk
import cleaning, xgb, ipdb, time

train_raw = pd.read_csv("train.csv")
test_raw  = pd.read_csv("test.csv")

train = cleaning.clean_data(train_raw)
test  = cleaning.clean_data(test_raw)

ipdb.set_trace()

# Identifying response, y, and features, X
y_train = train["y"]
x_train = train.drop(["ID","y"], axis = 1)
