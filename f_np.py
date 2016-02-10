# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 19:20:00 2016

@author: Jingtao
"""

import numpy as np
from sft import sft

def f_np(phi, lambd):
    """calculate f(lambda) based on the p order estimated AR models"""
    f = 1. / (np.absolute(1. - sft(phi, lambd)) ** 2)
    return f