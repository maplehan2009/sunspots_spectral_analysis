# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:39:32 2016

@author: Jingtao
"""

# import the packages
import numpy as np
import matplotlib.pyplot as plt
from autocovariance import autocovariance
from periodogram import periodogram
from f_np import f_np
from scipy import stats
from scipy.stats import probplot
from sft import sft

#------------------------------------------------------------------------------
# load data 
path = 'yearssn.dat'
raw_data = np.loadtxt(path)
sunspots = raw_data[:, 1]

# some basic statistical descriptions
size = sunspots.size
mean = np.mean(sunspots)
std = np.std(sunspots)
Q_0 = np.min(sunspots)
Q_25 = np.percentile(sunspots, 25)
Q_50 = np.percentile(sunspots, 50)
Q_75 = np.percentile(sunspots, 75)
Q_100 = np.max(sunspots)

print '\n' + '*' * 50
print 'Basic Statistical Description of the sunspots'
print '{0:10}    {1:10}'.format("attribute", "value")
print '{0:10}    {1:10f}'.format("count", size)
print '{0:10}    {1:10f}'.format("mean", mean)
print '{0:10}    {1:10f}'.format("std", std)
print '{0:10}    {1:10f}'.format("min", Q_0)
print '{0:10}    {1:10f}'.format("25%", Q_25)
print '{0:10}    {1:10f}'.format("50%", Q_50)
print '{0:10}    {1:10f}'.format("75%", Q_75)
print '{0:10}    {1:10f}'.format("max", Q_100)

# preprocessing of the data: centralize the data of sunspot
sunspots_centered = sunspots - mean

#------------------------------------------------------------------------------
# Question 4
# generate an array of autocovariance function, gamma[0], gamma[1]...gamma[p]
order_max = 10
phi_list = []
variance_list = []

for p in range(1, order_max+1):
    Gamma = [autocovariance(sunspots_centered, x) for x in range(p+1)]

    # solve the Yule-Walker equation with a certain order p
    vec_gamma = np.array(Gamma[1:])

    mat_gamma = np.zeros([p, p])
    for i in range(p):
        for j in range(p):
            h = np.abs(i - j)
            mat_gamma[i,j] = Gamma[h]

    coeff_phi = np.linalg.solve(mat_gamma, vec_gamma)
    variance_phi = Gamma[0] - np.dot(coeff_phi, vec_gamma)
    phi_list.append(coeff_phi)
    variance_list.append(variance_phi)
    
# print the results of phi and var    
for i in range(order_max):
    print '\n' + '*' * 50
    print "\nOrder p = " + str(i+1) + ":"
    print "Coefficients of prediction: "
    print phi_list[i]
    print "Error of prediction: " + str(variance_list[i])

#------------------------------------------------------------------------------
# Question 5
# calculate the periodigram and f_np using the SFT method

# Question 6
# test the similarity by plotting the QQPlot
# of I/f against the exponential distribution

# lambda is an interval of [0, PI], x_freq is the frequencey correspondant[0, 0.5]
lambd = np.linspace(0, np.pi, 200)
x_freq = lambd / (2 * np.pi)
# y_perio is the periodogram I of the centered sunspots number data
y_perio = periodogram(sunspots_centered, lambd)

print '\n' + '*' * 50
plt.semilogy(x_freq, y_perio)
plt.xlabel('frequency')
plt.ylabel('spectral density')
plt.title('Periodogram I')
plt.show()

# calculate the f_np with order p
for p in range(order_max):
    # phi is the prediction coefficients
    phi = phi_list[p]
    # sigma is the variance of the residues
    sigma = variance_list[p]
    # f is the estimated spectral density f_np
    f = f_np(phi, lambd)
    # lmd is the mean value of the exponential distribution associated
    lmd = sigma / (2 * np.pi)
    
    print '\n' + '*' * 50
    plt.semilogy(x_freq, f)
    plt.xlabel('frequency')
    plt.ylabel('estimated spectral density')
    plt.title('f_np P=' + str(p+1))
    plt.show()
    
    # calculate y_p which is the quotient of periodogram I and f_np
    # theoretically, it follows the exponential law    
    y_p = y_perio / f
    # QQPlot of I/f against an exponential distribution
    probplot(y_p, sparams=(0, lmd), dist='expon', plot=plt)
    plt.xlabel('Exp Distribution')
    plt.ylabel('I/f')
    plt.title('QQPlot between I/f and Exp Distribution P=' + str(p+1))
    plt.show()
    

#------------------------------------------------------------------------------
# Question 7
# train the AR model with the first half of the data
# From the previous test, we select p=9 as our order in this question.
p = 9;
Gamma = [autocovariance(sunspots_centered[:size/2], x) for x in range(p+1)]
vec_gamma = np.array(Gamma[1:])
mat_gamma = np.zeros([p, p])
for i in range(p):
    for j in range(p):
        h = np.abs(i - j)
        mat_gamma[i,j] = Gamma[h]

coeff_phi = np.linalg.solve(mat_gamma, vec_gamma)
variance_phi = Gamma[0] - np.dot(coeff_phi, vec_gamma)
    
# calculate the periodogram with the second half of the data
y_perio = periodogram(sunspots_centered[size/2:], lambd)

print '\n' + '*' *50
plt.semilogy(x_freq, y_perio)
plt.xlabel('frequency')
plt.ylabel('spectral density')
plt.title('Periodogram I (second half of the data)')
plt.show()

# calculate the f_np with the first half of the data
f = f_np(coeff_phi, lambd)

print '\n' + '*' * 50
plt.semilogy(x_freq, f)
plt.xlabel('frequency')
plt.ylabel('estimated spectral density')
plt.title('f_np P=9 (first half of the data)')
plt.show()

y_p = y_perio / f
lmd = variance_phi / (2 * np.pi)
probplot(y_p, sparams=(0, lmd), dist='expon', plot=plt)
plt.xlabel('Exp Distribution')
plt.ylabel('I/f')
plt.title('QQPlot between I/f and Exp Distribution P=9 with different training sets')
plt.show()