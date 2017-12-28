# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:23:13 2017

@author: mat.mathews
"""

import numpy as np
import matplotlib.pyplot as plt


filename = './init_field_hit.dat'
dataIn = np.loadtxt(filename,dtype=np.double)

N = 129
U   = np.empty((N,N,N),dtype=np.double)
V   = np.empty((N,N,N),dtype=np.double)
W   = np.empty((N,N,N),dtype=np.double)

for k in range(0,N):
    for j in range(0,N):
        for i in range(0,N):
            ii = k*N*N + j*N + i 
            U[i,j,k] = dataIn[ii,3]
            V[i,j,k] = dataIn[ii,4]
            W[i,j,k] = dataIn[ii,5]

#%%
uprime = 0
q = 0
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):       
            uprime += (U[i,j,k]**2 + V[i,j,k]**2 + W[i,j,k]**2)/3
            q      += (U[i,j,k]**2 + V[i,j,k]**2 + W[i,j,k]**2)     
            
uprime = uprime/(N-1)/(N-1)/(N-1)
q = q/(N-1)/(N-1)/(N-1)
uprime = np.sqrt(uprime)
q = np.sqrt(q)
q_old = q

k0 = 4
A = 16*np.sqrt(2/np.pi)*uprime**2/k0**5

#%%
uprimeGoal = np.sqrt(1/3)
Unew = U*uprimeGoal/uprime
Vnew = V*uprimeGoal/uprime
Wnew = W*uprimeGoal/uprime

#%%

uprime = 0
q = 0
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):       
            uprime += (Unew[i,j,k]**2 + Vnew[i,j,k]**2 + Wnew[i,j,k]**2)/3
            q      += (Unew[i,j,k]**2 + Vnew[i,j,k]**2 + Wnew[i,j,k]**2)     
            
uprime = uprime/(N-1)/(N-1)/(N-1)
q = q/(N-1)/(N-1)/(N-1)
uprime = np.sqrt(uprime)
q = np.sqrt(q)

#%%

k0 = 4
mu = 0.0069
dx = 2*np.pi/128
A = 16*np.sqrt(2/np.pi)*uprime**2/k0**5
kineng = (3*A/64)*np.sqrt(2*np.pi)*k0**5
Li = np.sqrt(2*np.pi)/k0
tau = np.sqrt(32/A)*(2*np.pi)**(1/4)*k0**(-7/2)
omega = (15*A/256)*np.sqrt(2*np.pi)*k0**7
kolEta = (2*omega*(mu)**-2)**(-1/4)
#%%
f = open('U_Mt0p3_N128_k4.dat','w');
g = open('V_Mt0p3_N128_k4.dat','w');
h = open('W_Mt0p3_N128_k4.dat','w');
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):
            f.write("".join([str(Unew[i,j,k]), "\n"]))
            g.write("".join([str(Vnew[i,j,k]), "\n"]))
            h.write("".join([str(Wnew[i,j,k]), "\n"]))
f.close()
g.close()
h.close()
