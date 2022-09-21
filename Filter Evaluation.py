# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 12:49:06 2022

@author: Mahmoud Saad
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import dot
from numpy import dot, sum, tile, linalg
from numpy.linalg import inv
from numpy.random import uniform
import scipy


from filterpy.monte_carlo import systematic_resample
from numpy.linalg import norm
from numpy.random import randn
import scipy.stats
from numpy.random import seed
from matplotlib.patches import Patch


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
x_ref= A*np.sin(omega*t+phase_shift) # Nonlinear state dynamics


# Measurment Model
gain=1
sigma=0.5
variance=sigma**2
offset=0
b= np.random.normal(0,sigma,size)
measurements=gain*x_ref+b+offset

#%% Kalman Filter

def kf_predict(X, P, A, Q, B, U):
 X = dot(A, X) + dot(B, U)
 P = dot(A, dot(P, A.T)) + Q
 return(X,P)

def kf_update(X, P, Y, H, R):
 IM = dot(H, X)
 IS = R + dot(H, dot(P, H.T))
 Z=dot(H.T, inv(IS))
 K = dot(P, Z)
 X = X + dot(K, (Y-IM))
 P = P - dot(K, dot(IS, K.T))
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

#%% Application of Kalman Filter 


dt = 0.1
# Initialization of state matrices
X = np.array([0.0])
P = np.diag([[10,10,10]])
A = np.array([[cos()]])
Q = np.array([[0]])
B = np.array([[0.1]])
U = np.array([1])

# Measurement matrices
Y = y[0]
H = np.array([[1]])
R = np.array([[10]])
# Number of iterations in Kalman Filter
N_iter = 1000

# Applying the Kalman Filter
Y_id=[]
for i in np.arange(0, N_iter):
  Y = y[i]
  (X, P) = kf_predict(X, P, A, Q, B, U)
  (X, P, K, IM, IS, LH) = kf_update(X, P, Y, H, R)
  Y_id.append(X[0])
  


#%% Particle Filter

def create_uniform_particles(y_range,N):
    particles = np.empty((N, 1))
    particles[:, 0] = uniform(y_range[0], y_range[1], size=N)
    return particles

def create_gaussian_particles(mean, std, N):
    particles =np.random.normal(0,2,1000)
   
    return particles

def predict(particles, u):

    # N = len(particles)
    # update particles position   
    particles += 0.1

def update(particles, weights, R, measurement):
      
 
   distance = np.linalg.norm(particles - measurement, axis=1)
   # print(len(distance))
   weights *= scipy.stats.norm(distance, R).pdf(distance)
   print(weights)

   weights += 1.e-300      # avoid round-off to zero
   weights /= sum(weights) # normalize

def estimate(particles, weights):
    """returns mean and variance of the weighted particles"""

    pos = particles
    mean = np.average(pos, weights=weights, axis=0)
    var  = np.average((pos - mean)**2, weights=weights, axis=0)
    return mean, var

def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights.resize(len(particles))
    weights.fill (1.0 / len(weights))
    
def neff(weights):
    return 1. / np.sum(np.square(weights))

#%% Application of the particle filter
N=1000
particles= np.random.normal(0,5,N)
weights= np.ones(N) /N
u=0
x_estimate=[]
k=0
for measurement in y: 
    
    predict(particles,u)
    u+=0.1
    # diff= particles - measurement
    # distance= np.linalg.norm(diff)
    weights *= scipy.stats.norm(measurement, 5).pdf(measurement) # How to update weight?
    weights += 1.e-300      # avoid round-off to zero
    weights /= sum(weights) # normalize weights
    
    # if neff(weights) < N/2:
    #         indexes = systematic_resample(weights)
    #         resample_from_index(particles, weights, indexes)
    #         assert np.allclose(weights, 1/N)
    
    mu, var = estimate(particles, weights)
    x_estimate.append(mu)
    

#%% Plotting yielded Signals

Measured_signal_color="red"
PF_curve_color="green"
KF_curve_color="black"


Measured_signal_Patch= Patch(color=Measured_signal_color, label='Disturbed  signal')
PF_Patch= Patch(color=PF_curve_color, label='Particle Filter')
KF_Patch= Patch(color=KF_curve_color, label='Kalman Filter')

patches=[Measured_signal_Patch,PF_Patch,KF_Patch]



plt.figure()
plt.plot(t,y,color=Measured_signal_color)
plt.plot(t,x_estimate, color=PF_curve_color)
plt.plot(t,Y_id,color=KF_curve_color)
plt.plot(t,x,color=KF_curve_color)

plt.title("Particle Filter vs Kalman Filter")
plt.grid(which='both')
plt.xlabel('time[s]')
plt.ylabel("R_on[ohm]")
plt.legend(handles=patches,loc='best')





#%% Bayes Filter 




