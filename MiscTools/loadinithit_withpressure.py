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
U   = np.empty((N-1,N-1,N-1),dtype=np.double)
V   = np.empty((N-1,N-1,N-1),dtype=np.double)
W   = np.empty((N-1,N-1,N-1),dtype=np.double)
dx = 2*np.pi/N
dy = 2*np.pi/N
dz = 2*np.pi/N

for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):
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
uprimeGoal = 1
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

U = Unew
V = Vnew
W = Wnew

#%%

#Need to get the source field for the poisson eqn.
#these are rough perturbations for the field, 2nd order central should be enough

S = np.empty((N-1,N-1,N-1),dtype=np.double)

Ux = np.empty((N-1,N-1,N-1),dtype=np.double)
Uy = np.empty((N-1,N-1,N-1),dtype=np.double)
Uz = np.empty((N-1,N-1,N-1),dtype=np.double)

Vx = np.empty((N-1,N-1,N-1),dtype=np.double)
Vy = np.empty((N-1,N-1,N-1),dtype=np.double)
Vz = np.empty((N-1,N-1,N-1),dtype=np.double)

Wx = np.empty((N-1,N-1,N-1),dtype=np.double)
Wy = np.empty((N-1,N-1,N-1),dtype=np.double)
Wz = np.empty((N-1,N-1,N-1),dtype=np.double)

for k in range(0,N-1):
    for j in range(0, N-1):
        for i in range(0, N-1):
            if i==0:
                Ux[0,j,k] = (U[1,j,k] - U[-1,j,k])/(2*dx)
                Vx[0,j,k] = (V[1,j,k] - V[-1,j,k])/(2*dx)
                Wx[0,j,k] = (W[1,j,k] - W[-1,j,k])/(2*dx)
            elif i==(N-2):
                Ux[-1,j,k] = (U[0,j,k] - U[-2,j,k])/(2*dx)
                Vx[-1,j,k] = (V[0,j,k] - V[-2,j,k])/(2*dx)
                Wx[-1,j,k] = (W[0,j,k] - W[-2,j,k])/(2*dx)
            else:
                Ux[i,j,k] = (U[i+1,j,k] - U[i-1,j,k])/(2*dx)
                Vx[i,j,k] = (V[i+1,j,k] - V[i-1,j,k])/(2*dx)
                Wx[i,j,k] = (W[i+1,j,k] - W[i-1,j,k])/(2*dx)    
                
            if j==0:
                Uy[i,0,k] = (U[i,1,k] - U[i,-1,k])/(2*dx)
                Vy[i,0,k] = (V[i,1,k] - V[i,-1,k])/(2*dx)
                Wy[i,0,k] = (W[i,1,k] - W[i,-1,k])/(2*dx)
            elif j==(N-2):
                Uy[i,-1,k] = (U[i,0,k] - U[i,-2,k])/(2*dx)
                Vy[i,-1,k] = (V[i,0,k] - V[i,-2,k])/(2*dx)
                Wy[i,-1,k] = (W[i,0,k] - W[i,-2,k])/(2*dx)
            else:
                Uy[i,j,k] = (U[i,j+1,k] - U[i,j-1,k])/(2*dx)
                Vy[i,j,k] = (V[i,j+1,k] - V[i,j-1,k])/(2*dx)
                Wy[i,j,k] = (W[i,j+1,k] - W[i,j-1,k])/(2*dx)   
                
            if k==0:
                Uz[i,j,0] = (U[i,j,1] - U[i,j,-1])/(2*dx)
                Vz[i,j,0] = (V[i,j,1] - V[i,j,-1])/(2*dx)
                Wz[i,j,0] = (W[i,j,1] - W[i,j,-1])/(2*dx)
            elif k==(N-2):
                Uz[i,j,-1] = (U[i,j,0] - U[i,j,-2])/(2*dx)
                Vz[i,j,-1] = (V[i,j,0] - V[i,j,-2])/(2*dx)
                Wz[i,j,-1] = (W[i,j,0] - W[i,j,-2])/(2*dx)
            else:
                Uz[i,j,k] = (U[i,j,k+1] - U[i,j,k-1])/(2*dx)
                Vz[i,j,k] = (V[i,j,k+1] - V[i,j,k-1])/(2*dx)
                Wz[i,j,k] = (W[i,j,k+1] - W[i,j,k-1])/(2*dx) 

S = -(Ux*Ux + Vy*Vy + Wz*Wz + 2*Uy*Vx + 2*Vz*Wy + 2*Uz*Wx)

#%% 

ptilde1 = np.empty((N-1,N-1,N-1),dtype=np.double)
ptilde2 = np.empty((N-1,N-1,N-1),dtype=np.double)
ppadtemp = np.empty((N+1,N+1,N+1),dtype=np.double)

r = np.empty((N-1,N-1,N-1),dtype=np.double)

ptilde1[:,:,:] = 0.0
ptilde2[:,:,:] = 0.0
r[:,:,:] = 0.0

omega = 2 / ( 1 + np.sin(np.pi/(N)) )

for kk in range(0,2000):

    ppadtemp[:,:,:] = 0.0
    ppadtemp[1:-1,1:-1, 1:-1]   = ptilde1
    ppadtemp[0   ,1:-1, 1:-1]   = ptilde1[-1,:,:]
    ppadtemp[-1  ,1:-1, 1:-1]   = ptilde1[ 0,:,:]
    ppadtemp[1:-1,0   , 1:-1]   = ptilde1[:,-1,:]
    ppadtemp[1:-1,-1  , 1:-1]   = ptilde1[:, 0,:]
    ppadtemp[1:-1,1:-1,    0]   = ptilde1[:,:,-1]
    ppadtemp[1:-1,1:-1,   -1]   = ptilde1[:,:, 0]

    
    r[:,:,:] = 0.0
    r -= S*dx*dx
    #r[i,j,k]   -= 0 #6*ptilde1[i,j,k]
    r[:,:,:] += ppadtemp[0:-2,1:-1,1:-1]
    r[:,:,:] += ppadtemp[2:  ,1:-1,1:-1]
    
    r[:,:,:] += ppadtemp[1:-1,0:-2,1:-1]
    r[:,:,:] += ppadtemp[1:-1,2:  ,1:-1]    
    
    r[:,:,:] += ppadtemp[1:-1,1:-1,0:-2]
    r[:,:,:] += ppadtemp[1:-1,1:-1,2:]

                       
    ptilde2 = (1/6)*r
    res_norm = np.sum((ptilde2-ptilde1)**2)
    ptilde1 = ptilde2
    
    print(kk)
    print(res_norm)
    
    if res_norm < 0.001:
        break
    
    
#%%

Pxx = np.empty((N-1,N-1,N-1),dtype=np.double)
Pyy = np.empty((N-1,N-1,N-1),dtype=np.double)
Pzz = np.empty((N-1,N-1,N-1),dtype=np.double)
PS  = np.empty((N-1,N-1,N-1),dtype=np.double)


for k in range(0,N-1):
    for j in range(0, N-1):
        for i in range(0, N-1):
            if i==0:
                Pxx[0,j,k] = (ptilde1[1,j,k] -2*ptilde1[0,j,k] + ptilde1[-1,j,k])/(dx*dx)
            elif i==(N-2):
                Pxx[-1,j,k] = (ptilde1[0,j,k] -2*ptilde1[-1,j,k] + ptilde1[-2,j,k])/(dx*dx)
            else:
                Pxx[i,j,k] = (ptilde1[i+1,j,k] -2*ptilde1[i,j,k] + ptilde1[i-1,j,k])/(dx*dx)
                
            if j==0:
                Pyy[i,0,k] = (ptilde1[i,1,k] -2*ptilde1[i,0,k] + ptilde1[i,-1,k])/(dx*dx)
            elif j==(N-2):
                Pyy[i,-1,k] = (ptilde1[i,0,k] -2*ptilde1[i,-1,k] + ptilde1[i,-2,k])/(dx*dx)
            else:
                Pyy[i,j,k] = (ptilde1[i,j+1,k] -2*ptilde1[i,j,k] + ptilde1[i,j-1,k])/(dx*dx)
                
            if k==0:
                Pzz[i,j,0] = (ptilde1[i,j,1] -2*ptilde1[i,j,0] + ptilde1[i,j,-1])/(dx*dx)
            elif k==(N-2):
                Pzz[i,j,-1] = (ptilde1[i,j,0] -2*ptilde1[i,j,-1] + ptilde1[i,j,-2])/(dx*dx)
            else:
                Pzz[i,j,k] = (ptilde1[i,j,k+1] -2*ptilde1[i,j,k] + ptilde1[i,j,k-1])/(dx*dx)

PS = Pxx + Pyy + Pzz

#%%

#Chop the domain into three chunks
#Chunk1
U1 = Unew[:,:,0:42]
V1 = Vnew[:,:,0:42]
W1 = Wnew[:,:,0:42]
P1 = ptilde1[:,:,0:42]

#Chunk2
U2 = Unew[:,:,42:84]
V2 = Vnew[:,:,42:84]
W2 = Wnew[:,:,42:84]
P2 = ptilde1[:,:,42:84]

#Chunk3
U3 = Unew[:,:,84:126]
V3 = Vnew[:,:,84:126]
W3 = Wnew[:,:,84:126]
P3 = ptilde1[:,:,84:126]

totalX = 256
currentX = 128*3
totalOverlap = currentX - totalX

Ufinal = np.empty((totalX,N-1,42),dtype=np.double)
Vfinal = np.empty((totalX,N-1,42),dtype=np.double)
Wfinal = np.empty((totalX,N-1,42),dtype=np.double)
Pfinal = np.empty((totalX,N-1,42),dtype=np.double)

Ufinal[42:86,:,:] = U1[42:86,:,:]
Vfinal[42:86,:,:] = V1[42:86,:,:]
Wfinal[42:86,:,:] = W1[42:86,:,:]
Pfinal[42:86,:,:] = P1[42:86,:,:]

Ufinal[128:172,:,:] = U2[42:86,:,:]
Vfinal[128:172,:,:] = V2[42:86,:,:]
Wfinal[128:172,:,:] = W2[42:86,:,:]
Pfinal[128:172,:,:] = P2[42:86,:,:]

Ufinal[214:256,:,:] = U3[42:84,:,:]
Vfinal[214:256,:,:] = V3[42:84,:,:]
Wfinal[214:256,:,:] = W3[42:84,:,:]
Pfinal[214:256,:,:] = P3[42:84,:,:]
for i in range(0,42):
    theta = (np.pi/2.0)*float(i)/41.0
    Ufinal[i,:,:] = np.cos(theta)*U3[84+i,:,:] + np.sin(theta)*U1[i,:,:]
    Vfinal[i,:,:] = np.cos(theta)*V3[84+i,:,:] + np.sin(theta)*V1[i,:,:]
    Wfinal[i,:,:] = np.cos(theta)*W3[84+i,:,:] + np.sin(theta)*W1[i,:,:]
    Pfinal[i,:,:] = np.cos(theta)*P3[84+i,:,:] + np.sin(theta)*P1[i,:,:]
    
    Ufinal[86+i,:,:] = np.cos(theta)*U1[86+i,:,:] + np.sin(theta)*U2[i,:,:]
    Vfinal[86+i,:,:] = np.cos(theta)*V1[86+i,:,:] + np.sin(theta)*V2[i,:,:]
    Wfinal[86+i,:,:] = np.cos(theta)*W1[86+i,:,:] + np.sin(theta)*W2[i,:,:]
    Pfinal[86+i,:,:] = np.cos(theta)*P1[86+i,:,:] + np.sin(theta)*P2[i,:,:]    

    Ufinal[172+i,:,:] = np.cos(theta)*U2[86+i,:,:] + np.sin(theta)*U3[i,:,:]
    Vfinal[172+i,:,:] = np.cos(theta)*V2[86+i,:,:] + np.sin(theta)*V3[i,:,:]
    Wfinal[172+i,:,:] = np.cos(theta)*W2[86+i,:,:] + np.sin(theta)*W3[i,:,:]
    Pfinal[172+i,:,:] = np.cos(theta)*P2[86+i,:,:] + np.sin(theta)*P3[i,:,:]    

#%%
f = open('U_Mt0p3_N128_k8.dat','w');
g = open('V_Mt0p3_N128_k8.dat','w');
h = open('W_Mt0p3_N128_k8.dat','w');
pp = open('Pdiv_gamma_rho_M2_N128_k8.dat','w')
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):
            f.write("".join([str(Unew[i,j,k]), "\n"]))
            g.write("".join([str(Vnew[i,j,k]), "\n"]))
            h.write("".join([str(Wnew[i,j,k]), "\n"]))
            pp.write("".join([str(ptilde1[i,j,k]),"\n]))
f.close()
g.close()
h.close()

#%%

f = open('U_uprime1_N128_k8_256x128x42.dat','w');
g = open('V_uprime1_N128_k8_256x128x42.dat','w');
h = open('W_uprime1_N128_k8_256x128x42.dat','w');
pp = open('P_uprime1_N128_k8_256x128x42.dat','w')
for k in range(0,42):
    for j in range(0,128):
        for i in range(0,256):
            f.write("".join([str(Ufinal[i,j,k]), "\n"]))
            g.write("".join([str(Vfinal[i,j,k]), "\n"]))
            h.write("".join([str(Wfinal[i,j,k]), "\n"]))
            pp.write("".join([str(Pfinal[i,j,k]), "\n"]))
f.close()
g.close()
h.close()
pp.close()