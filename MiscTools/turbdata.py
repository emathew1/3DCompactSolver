#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 11:56:17 2017

@author: edwin
"""

import numpy as np
import matplotlib.pyplot as plt

filename = './turbdata.out'
dataIn = np.loadtxt(filename,dtype=np.double)


time = dataIn[:,0]
tau  = dataIn[:,1]
taylorReyn = dataIn[:,2]
turbMach   = dataIn[:,3]
meanTurbDiss = dataIn[:,4]
uprime = dataIn[:,5]
kolNu  = dataIn[:,6]
rhoprime = dataIn[:,7]
dilprime = dataIn[:,8]
vortprime = dataIn[:,9]
meanKineticEng = dataIn[:,10]

#%%
filename = './D4_Mt.dat'
dataIn = np.loadtxt(filename,dtype=np.double)
Mt_time = dataIn[:,0]
Mt_Mt   = dataIn[:,1]

filename = './D4_K.dat'
dataIn = np.loadtxt(filename,dtype=np.double)
K_time = dataIn[:,0]
K_K   = dataIn[:,1]

filename = './D4_Rhoprime.dat'
dataIn = np.loadtxt(filename,dtype=np.double)
Rho_time = dataIn[:,0]
Rho_Rho   = dataIn[:,1]



#%%
plt.plot(time/tau[0],dilprime/vortprime[0])

#%%
plt.plot(Rho_time, Rho_Rho, 'k-', time/tau[0],rhoprime/turbMach[0]**2, 'b-')

#%%
plt.plot(K_time, K_K,'k-',time/tau[0],meanKineticEng/meanKineticEng[0],'b-')

#%%
plt.loglog(time/tau[0],taylorReyn)

#%%
plt.plot(time/tau[0],dilprime)

#%%
plt.plot(Mt_time, Mt_Mt,'k-o',time/tau[0],turbMach/turbMach[0],'b-')

#%%
filename = './turbdata_D4k4.out'
dataIn = np.loadtxt(filename,dtype=np.double)


time = dataIn[:,0]
tau  = dataIn[:,1]
taylorReyn = dataIn[:,2]
turbMach   = dataIn[:,3]
meanTurbDiss = dataIn[:,4]
uprime = dataIn[:,5]
kolNu  = dataIn[:,6]
rhoprime = dataIn[:,7]
dilprime = dataIn[:,8]
vortprime = dataIn[:,9]
meanKineticEng = dataIn[:,10]

filename = './turbdata_D4k8.out'
dataIn = np.loadtxt(filename,dtype=np.double)

time2 = dataIn[:,0]
tau2  = dataIn[:,1]
taylorReyn2 = dataIn[:,2]
turbMach2   = dataIn[:,3]
meanTurbDiss2 = dataIn[:,4]
uprime2 = dataIn[:,5]
kolNu2  = dataIn[:,6]
rhoprime2 = dataIn[:,7]
dilprime2 = dataIn[:,8]
vortprime2 = dataIn[:,9]
meanKineticEng2 = dataIn[:,10]