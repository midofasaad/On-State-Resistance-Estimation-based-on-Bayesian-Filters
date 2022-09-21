# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:58:13 2022

@author: Mahmoud

Description: Variants of the Particle Filter are applied to determine different best estimates for a reference signal with known system dynamics.
             The reference signal reperesents the ideal signal, with which the estimates from the filters are compared.
             Measurments are simulated in time interval [0,tmax] through scaling of the state with factor a and the addition of 
             gaussian noise. Finally, a Graphical comparsion through Matplotlib is  preformed to assess the accuracy of each variant 
             in light of the reference signal.
"""



import numpy as np
import matplotlib.pyplot as plt
import scipy
from filterpy.monte_carlo import systematic_resample
import scipy.stats
from scipy.interpolate import interp1d,BarycentricInterpolator,InterpolatedUnivariateSpline
import scipy.optimize as optimization




#%% Particle Filter-Based Estimation
def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights.resize(len(particles))
    weights.fill (1.0 / len(weights))

def estimate(particles, weights):
    """returns mean and variance of the weighted particles"""

    pos = particles
    mean = np.average(pos, weights=weights, axis=0)
    var  = np.average((pos - mean)**2, weights=weights, axis=0)
    return mean, var

def neff(weights):
    return 1. / np.sum(np.square(weights))

def determine_parameters(t,x): 
    parameters=[0,0,0]
    if len(t)>=3: 
       
        X=np.array([np.ones_like(t),t,t**2],dtype='float').T
        # # F= X.T@X
        # # F=np.array(F,dtype='float')
        F=np.linalg.pinv(X)
        # F_inv=np.linalg.inv(F)
        # E=np.matmul(F_inv,X.T)
        parameters=F@x
        
        return parameters
    else: 
        return parameters
    
#%% Smoothing algorithms   
def moving_average_algorithm(n,t_slice,measurements,plot):
    size_y=len(measurements)
    y=np.empty(size_y)
    k=0
    factor=2*(n)+1
    for x in measurements:
        if k<size_y-1:
            if k>=n and size_y-1-k>=n:
                y[k]=sum(measurements[k-n:k+(n+1)])/factor
            else: 
                y[k]=x
        k+=1
    if plot==True:    
        plt.figure()
        plt.plot(t_slice,y,color='red')    
        plt.plot(t_slice,measurements)
    
    return y


#%%  Estimation method
def gaus_estimate(A_estimate,phase_estimate,c):
     c=lamda*c+0.5
     e=y-A_estimate*np.sin((omega/fs)*k+phase_estimate)
     phase_estimate_prev=phase_estimate
     phase_estimate=phase_estimate+np.cos((omega/fs)*k+phase_estimate)*(e/(A_estimate*c))
     A_estimate=A_estimate+np.sin((omega/fs)*k+phase_estimate_prev)*(e/c)
   
    





#%% Experimental Area



# Process  Model
delta_t=100e-6
size=1000
tmax=size*delta_t
freq=100
omega=2*np.pi*freq
A=10 #Amplitude
phase_shift=0
t=np.linspace(0,tmax,num=size)
x_ref= A*np.cos(omega*t+phase_shift) # Nonlinear state dynamics


# Measurment Model
gain=1
sigma=0.5
variance=sigma**2
offset=10
b= np.random.normal(0,sigma,size)
measurements=gain*x_ref+b+offset


# Parameters of the Particle Filter
a=-20
b=20
N=500
part=np.random.uniform(a,b,N)
part_2=np.random.uniform(a,b,N)

weight=np.ones(N)/N 
weight_2=np.ones(N)/N
 
z=delta_t*omega
r=1/((z)**2+1)

x_estimate=np.empty((size,1))
x_estimate_2=np.empty((size,1))

x_variance=np.empty((size,1))
x_variance_2=np.empty((size,1))


mu, var = estimate(part, weight)
x_estimate[0]=0
x_estimate_2[0]=0
# particles=np.empty((size,N))
# particles[:,0]=part
fitting_points=[]
k=1 
part_k=np.ones_like(part)*0 
variance_estimates=[]
sigma_estimates=[]
offset_estimate=0
offset_estimates=[]
online_measurements=[]
mu=0

# Amplitude and phase shift intialization factors
A_estimate=1
A_estimates=[]
phase_estimate=2
phase_estimates=[]
gaus_estimates=np.empty((size,1))
c=1
lamda=0.999
fs=1/delta_t
i=0

# Predict next state based on state dynamics
for y in measurements:
     online_measurements.append(y)
     offset_estimate=(y+(k-1)*offset_estimate)/k
     diff_square=sum((online_measurements-offset_estimate-mu)**2)
     variance_estimate=(1/len(online_measurements))*diff_square
     sigma_estimate=np.sqrt(variance_estimate)
     
     variance_estimates.append(variance_estimate)
     sigma_estimates.append(sigma_estimate)
     offset_estimates.append(offset_estimate)
     
     delta_offset=abs(offset_estimate-offset_estimates[k-2])
     y=y-offset_estimate
     
     c=lamda*c+0.5
     e=y-A_estimate*np.sin((omega/fs)*k+phase_estimate)
     phase_estimate_prev=phase_estimate
     phase_estimate=phase_estimate+np.cos((omega/fs)*k+phase_estimate)*(e/(A_estimate*c))
     A_estimate=A_estimate+np.sin((omega/fs)*k+phase_estimate_prev)*(e/c)
     
     A_estimates.append(A_estimate)
     phase_estimates.append(phase_estimate)
     
     gaus_estimate=A_estimate*np.sin((omega/fs)*k+phase_estimate)
     gaus_estimates[k]=gaus_estimate
     
     
           
     # if t[k]>0.02 and i==0:
     #     part=np.random.uniform(-2*A_estimate,2*A_estimate,N)
     #     weight=np.ones(N)/N 
     #     i=i+1
       
     
     # part=part+(t[k]-t[k-1])*(omega)*(A_estimate)*np.cos(omega*t[k-1])
     
    
     # part=r*(part+abs(A_estimate)*omega*np.cos(omega*t[k])*delta_t)
     # part_2=r*(part_2-z*part)
     
     part=np.cos(omega/fs)*part+(np.sin(omega/fs)/omega)*part_2
     part_2=-np.sin(omega/fs)*omega*part+np.cos(omega/fs)*part_2
     
     
     # particles[:,k]=part
     
     
     diff=(-y+part)/(sigma)
    
     weight*=scipy.stats.norm.pdf(diff)/(sigma)
     weight += 1.e-300
     weight/=sum(weight)
     
     N_eff=neff(weight)
    
     if N_eff<0.8*N: 
         indexes = systematic_resample(weight)
         resample_from_index(part, weight, indexes)
         assert np.allclose(weight, 1/N)
         
   
     mu, var = estimate(part, weight)
     
     
     k=k+1
     if k>=1000:
               break
     
         
     
     x_estimate[k]=mu
     x_variance[k]=var
     
     # x_estimate_2[k]=mu_2
     # x_variance_2[k]=var_2
     
   


    
    



#%% Plotting
plt.figure()
plt.plot(t,measurements,color='blue',label='Measurment')
#plt.plot(t,x_estimate,color='black',label='PF Estimate')
#plt.plot(t,x_estimate_2,color='yellow',label='PF Estimate_cos')
plt.plot(t,x_ref,color='red',label='desired signal')
#plt.plot(t,gaus_estimates,color="green",label='gaus-recursive estimate')
plt.xlabel('time[ms]')
plt.ylabel('Voltage[v]')
plt.legend()
plt.grid()

plt.figure()
plt.plot(t,x_variance)
plt.xlabel('time[ms]')
plt.ylabel('Absolute Error[v]')
plt.grid()
     
#plt.figure()
#plt.plot(t[:-1],parameters_mega)
 



     
     
     




