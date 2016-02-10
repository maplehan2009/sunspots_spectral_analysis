# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 15:07:07 2016

@author: Jingtao
"""

import numpy as np

def autocovariance(array, h):
    """array is the sample series
    h is the input variable of the autocovariance function
    Here, for simplicity, we assume the element of array, X_n has one dimension
    And the value of X_n is real number"""
    
    n = array.size
    u = np.mean(array)
    if (h > n-1) | (h < -(n-1)):
        return 0
    else:
        h = np.abs(h)
        vec_1 = array[h:] - u
        vec_2 = array[:n-h] - u
        return np.dot(vec_1, vec_2) / n   