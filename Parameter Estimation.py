# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:58:44 2022

@author: Mahmoud Saad

Core-Function: Visualize a normal Distribution with arbitrary mean and variance, such that
             parameter estimation methods (e.g maximum liklihood) are applied to determine
             the parameters of the distribution.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats


variance=0.5 
sigma=np.sqrt(variance)
mu= 10 
size=1000

interval=np.linspace(0,20,size)
prob_gauss=scipy.stats.norm.pdf(interval-mu/sigma)/sigma 

N=500
summ=0
online_measurements=[]
mu_estimates=[]
variance_estimates=[]
sigma_estimates=[]

for i in range(0,N+1): 
    measurement=sigma*np.random.randn(1)+mu
    online_measurements.append(measurement)
    mu_estimate=1/len(online_measurements)*sum(online_measurements)
    mu_estimates.append(mu_estimate)
    
    diff_square=sum((online_measurements-np.array(mu))**2)[0]
    variance_estimate=1/len(online_measurements)*diff_square
    sigma_estimate=np.sqrt(variance_estimate)
    
    variance_estimates.append(variance_estimate)
    sigma_estimates.append(sigma_estimate)
    
    
    
    
    
prob_gauss_estimate=scipy.stats.norm.pdf(interval-mu_estimate/sigma)/sigma

plt.figure()
plt.plot(interval,prob_gauss,label='Original Gaussian')
plt.plot(interval,prob_gauss_estimate,label='Estimated Gauss')
plt.xlabel('Possible Values of Random Variable X')
plt.ylabel('Associated Probability')
plt.legend()
plt.grid()


plt.figure()
plt.plot(range(0,N+1),mu_estimates,label='Estimate mean')
plt.plot(range(0,N+1),mu*np.ones_like(mu_estimates),label='Original mean')
plt.xlabel('Number of Samples [N]')
plt.ylabel('Mean')
plt.legend()
plt.grid()


plt.figure()
plt.plot(range(0,N+1),variance_estimates,label='Estimated Variance')
plt.plot(range(0,N+1),variance*np.ones_like(variance_estimates),label='Original Varience')
plt.xlabel('Number of Samples [N]')
plt.ylabel('Variance')
plt.legend()
plt.grid()







