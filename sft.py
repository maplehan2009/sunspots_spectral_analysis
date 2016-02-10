# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 18:48:16 2016

@author: Jingtao
"""
import numpy as np

def sft(X, lambd):
    """Slow Fourier Transform Function
    I do not use the FFT because of the output number.
    e.g. When the input has just one number, FFT will give the output
    with also only one number. I want the lambd spreaded from 0 to PI.
    So I write the Fourier Transform myself"""
    n = X.size
    ans = np.zeros(lambd.size)
    for i in range(n):
        ans = ans + (X[i] * np.cos((i+1) * lambd) - 1j * X[i] * np.sin((i+1) * lambd))
    return ans