# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 18:36:58 2016

@author: Jingtao
"""

import numpy as np
from sft import sft
  
def periodogram(X, lambd):
    """calculate the periodogram of a time series"""
    n = X.size
    I = np.absolute(sft(X, lambd)) ** 2 / (2 * np.pi * n)
    I[0] = 0
    return I    