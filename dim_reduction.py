import pandas as pd
import numpy as np
from sklearn.decomposition import PCA, FastICA

def _pca_reduction(X):
    pca = PCA(n_components = n_components)
    return pca.fit_transform(X)

def _ica_reduction(X, n_components = 10):
    ica = FastICA(n_components = n_components)
    return ica.fit_transform(X)

def reduce(X, reduction_type, n_componenets = 10):
    # Initialising output variable
    output = None
    reduction_type = reduction_type.upper()
    
    # Switch between types of reduction
    if (reduction_type == "PCA"):
        output = _pca_reduction(X, y, n_components)
    elif (reduction_type == "ICA"):
        output = _ica_reduction(X, y, n_components)
    else:
        raise ValueError("Reduction type %s has not been implemented yet" % reduction_type)
    
    return output
