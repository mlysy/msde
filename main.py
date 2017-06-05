import numpy as np
import pandas as pd
import sklearn as sk
import cleaning, xgb, ipdb, time

train_data = pd.read_csv("train.csv")
test_data  = pd.read_csv("test.csv")

train_clean = cleaning.clean_data(train_data)
test_clean  = cleaning.clean_data(test_data)

ipdb.set_trace()
