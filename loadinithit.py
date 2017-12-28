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

Unew = U*(0.3/q)
Vnew = V*(0.3/q)
Wnew = V*(0.3/q)

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
f = open('U_Mt0p3_N128_k8.dat','w');
g = open('V_Mt0p3_N128_k8.dat','w');
h = open('W_Mt0p3_N128_k8.dat','w');
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):
            f.write("".join([str(Unew[i,j,k]), "\n"]))
            g.write("".join([str(Vnew[i,j,k]), "\n"]))
            h.write("".join([str(Wnew[i,j,k]), "\n"]))
f.close()
g.close()
h.close()
