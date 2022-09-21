# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 14:35:13 2022

@author: Mahmoud Saad
Description: Application of the Kalman Filter for the extraction of a
             sin signal from gaussian noised signal
"""


from numpy.random import randn
from filterpy.kalman import UnscentedKalmanFilter
from filterpy.common import Q_discrete_white_noise
from filterpy.kalman import JulierSigmaPoints
import numpy as np
import matplotlib.pyplot as plt
import scipy
from filterpy.monte_carlo import systematic_resample
import scipy.stats
from scipy.interpolate import interp1d,BarycentricInterpolator,InterpolatedUnivariateSpline
import scipy.optimize as optimization

#Linear Kalman Filter
from numpy import dot
from numpy import dot, sum, tile, linalg
from numpy.linalg import inv
from numpy.random import uniform
import scipy




#%%  Linear Kalman Filter

def kf_predict(X, P, A, Q, B, U):
 X = A@X + B*U
 P = A@P@A.T + Q
 return(X,P)

def kf_update(X, P, Y, H, R):
 IM = H@X
 IS = R + H@P@H.T
 Z=H@np.linalg.pinv(IS)
 K = P@Z.T
 X = X + K@(Y-IM).T
 print(X)
 P = P - K.T@IS@K
 LH = np.array([[1]])
 return (X,P,K,IM,IS,LH)

def gauss_pdf(X, M, S):
    if M.shape()[1] == 1:
         DX = X - tile(M, X.shape()[1])
         E = 0.5 * sum(DX * (dot(inv(S), DX)), axis=0)
         E = E + 0.5 * M.shape()[0] * np.log(2 * np.pi) + 0.5 * np.log(np.det(S))
         P = np.exp(-E)
    elif X.shape()[1] == 1:
         DX = tile(X, M.shape()[1])- M
         E = 0.5 * sum(DX * (dot(inv(S), DX)), axis=0)
         E = E + 0.5 * M.shape()[0] * np.log(2 * np.pi) + 0.5 * np.log(np.det(S))
         P = np.exp(-E)
    else:
         DX = X-M
         E = 0.5 * dot(DX.T, dot(inv(S), DX))
         E = E + 0.5 * M.shape()[0] * np.log(2 * np.pi) + 0.5 * np.log(np.det(S))
         P = np.exp(-E)
    return (P[0],E[0]) 

#%% Process and Measurement Model 

# Reference Signal
delta_t=100e-6
size=1000
tmax=size*delta_t
freq=100
omega=2*np.pi*freq
A=5 #Amplitude
phase_shift=0
t=np.linspace(0,tmax,num=size)
x_ref= A*np.cos(omega*t+phase_shift) # Nonlinear state dynamics


# Measurment Model
gain=1
sigma=0.5
variance=sigma**2
offset=2
b= np.random.normal(0,sigma,size)
measurements=gain*x_ref+b+offset



#%% Application of Kalman Filter 
dt = 0.1
fs=size/tmax
omega_k=omega/fs
# Initialization of state matrices
X = np.array([[0],[0]])
P = np.array([[10000,0],
              [0,10000]])
A = np.array([[np.cos(omega), -np.sin(omega)],
              [np.sin(omega), np.cos(omega)]])
Q = np.array([[5000,0],
              [0,5000]])
B = np.array([[0],[0]])
U = 0

# Measurement matrices
H = np.array([[1,0]])
R = np.array([[0,0],
              [0,0]])
# Number of iterations in Kalman Filter
N_iter = 1000

# Applying the Kalman Filter
Kf_estimate=[]
offset_estimate=[]
for y in measurements:
  Y= np.array([[y]])
  (X, P) = kf_predict(X, P, A, Q, B, U)
  (X, P, K, IM, IS, LH) = kf_update(X, P, Y, H, R)
  Kf_estimate.append(X[0])
  offset_estimate.append(X[1])

#%% Plotting
plt.figure()
plt.plot(t,measurements,color='blue',label='Measurment')
plt.plot(t,Kf_estimate,color='black',label='KF Estimate')
plt.plot(t,x_ref,color='red',label='real signal')
plt.xlabel('time[ms]')
plt.ylabel('Voltage[v]')
plt.legend()
plt.grid()


plt.figure()
plt.plot(t,offset_estimate,color='black',label='Offset Estimate')
plt.plot(t,offset*np.ones_like(t),color='red',label='real signal')
plt.xlabel('time[ms]')
plt.ylabel('Voltage[v]')
plt.legend()
plt.grid()