# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:23:13 2017

@author: mat.mathews
"""

import numpy as np
import matplotlib.pyplot as plt


filename1 = './init_field_hit.dat'
filename2 = './init_field_hit_2.dat'
dataIn1 = np.loadtxt(filename1,dtype=np.double)
dataIn2 = np.loadtxt(filename2,dtype=np.double)

N = 129
U1   = np.empty((N-1,N-1,N-1),dtype=np.double)
V1   = np.empty((N-1,N-1,N-1),dtype=np.double)
W1   = np.empty((N-1,N-1,N-1),dtype=np.double)
U2   = np.empty((N-1,N-1,N-1),dtype=np.double)
V2   = np.empty((N-1,N-1,N-1),dtype=np.double)
W2   = np.empty((N-1,N-1,N-1),dtype=np.double)
dx = 2*np.pi/N
dy = 2*np.pi/N
dz = 2*np.pi/N

for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):
            ii = k*N*N + j*N + i 
            U1[i,j,k] = dataIn1[ii,3]
            V1[i,j,k] = dataIn1[ii,4]
            W1[i,j,k] = dataIn1[ii,5]
            U2[i,j,k] = dataIn2[ii,3]
            V2[i,j,k] = dataIn2[ii,4]
            W2[i,j,k] = dataIn2[ii,5]            

#%%
uprime1 = 0
uprime2 = 0
q1 = 0
q2 = 0
for k in range(0,N-1):
    for j in range(0,N-1):
        for i in range(0,N-1):       
            uprime1 += (U1[i,j,k]**2 + V1[i,j,k]**2 + W1[i,j,k]**2)/3
            q1      += (U1[i,j,k]**2 + V1[i,j,k]**2 + W1[i,j,k]**2)     
            uprime2 += (U2[i,j,k]**2 + V2[i,j,k]**2 + W2[i,j,k]**2)/3
            q2      += (U2[i,j,k]**2 + V2[i,j,k]**2 + W2[i,j,k]**2)  
            
uprime1 = uprime1/(N-1)/(N-1)/(N-1)
q1 = q1/(N-1)/(N-1)/(N-1)
uprime1 = np.sqrt(uprime1)
q1 = np.sqrt(q1)

uprime2 = uprime2/(N-1)/(N-1)/(N-1)
q2 = q2/(N-1)/(N-1)/(N-1)
uprime2 = np.sqrt(uprime2)
q2 = np.sqrt(q2)

#%%
uprimeGoal = 1
U1 = U1*uprimeGoal/uprime1
V1 = V1*uprimeGoal/uprime1
W1 = W1*uprimeGoal/uprime1

U2 = U2*uprimeGoal/uprime2
V2 = V2*uprimeGoal/uprime2
W2 = W2*uprimeGoal/uprime2


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
                Ux[0,j,k] = (U1[1,j,k] - U1[-1,j,k])/(2*dx)
                Vx[0,j,k] = (V1[1,j,k] - V1[-1,j,k])/(2*dx)
                Wx[0,j,k] = (W1[1,j,k] - W1[-1,j,k])/(2*dx)
            elif i==(N-2):
                Ux[-1,j,k] = (U1[0,j,k] - U1[-2,j,k])/(2*dx)
                Vx[-1,j,k] = (V1[0,j,k] - V1[-2,j,k])/(2*dx)
                Wx[-1,j,k] = (W1[0,j,k] - W1[-2,j,k])/(2*dx)
            else:
                Ux[i,j,k] = (U1[i+1,j,k] - U1[i-1,j,k])/(2*dx)
                Vx[i,j,k] = (V1[i+1,j,k] - V1[i-1,j,k])/(2*dx)
                Wx[i,j,k] = (W1[i+1,j,k] - W1[i-1,j,k])/(2*dx)    
                
            if j==0:
                Uy[i,0,k] = (U1[i,1,k] - U1[i,-1,k])/(2*dx)
                Vy[i,0,k] = (V1[i,1,k] - V1[i,-1,k])/(2*dx)
                Wy[i,0,k] = (W1[i,1,k] - W1[i,-1,k])/(2*dx)
            elif j==(N-2):
                Uy[i,-1,k] = (U1[i,0,k] - U1[i,-2,k])/(2*dx)
                Vy[i,-1,k] = (V1[i,0,k] - V1[i,-2,k])/(2*dx)
                Wy[i,-1,k] = (W1[i,0,k] - W1[i,-2,k])/(2*dx)
            else:
                Uy[i,j,k] = (U1[i,j+1,k] - U1[i,j-1,k])/(2*dx)
                Vy[i,j,k] = (V1[i,j+1,k] - V1[i,j-1,k])/(2*dx)
                Wy[i,j,k] = (W1[i,j+1,k] - W1[i,j-1,k])/(2*dx)   
                
            if k==0:
                Uz[i,j,0] = (U1[i,j,1] - U1[i,j,-1])/(2*dx)
                Vz[i,j,0] = (V1[i,j,1] - V1[i,j,-1])/(2*dx)
                Wz[i,j,0] = (W1[i,j,1] - W1[i,j,-1])/(2*dx)
            elif k==(N-2):
                Uz[i,j,-1] = (U1[i,j,0] - U1[i,j,-2])/(2*dx)
                Vz[i,j,-1] = (V1[i,j,0] - V1[i,j,-2])/(2*dx)
                Wz[i,j,-1] = (W1[i,j,0] - W1[i,j,-2])/(2*dx)
            else:
                Uz[i,j,k] = (U1[i,j,k+1] - U1[i,j,k-1])/(2*dx)
                Vz[i,j,k] = (V1[i,j,k+1] - V1[i,j,k-1])/(2*dx)
                Wz[i,j,k] = (W1[i,j,k+1] - W1[i,j,k-1])/(2*dx) 

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


p1 = ptilde1

#%%

for k in range(0,N-1):
    for j in range(0, N-1):
        for i in range(0, N-1):
            if i==0:
                Ux[0,j,k] = (U2[1,j,k] - U2[-1,j,k])/(2*dx)
                Vx[0,j,k] = (V2[1,j,k] - V2[-1,j,k])/(2*dx)
                Wx[0,j,k] = (W2[1,j,k] - W2[-1,j,k])/(2*dx)
            elif i==(N-2):
                Ux[-1,j,k] = (U2[0,j,k] - U2[-2,j,k])/(2*dx)
                Vx[-1,j,k] = (V2[0,j,k] - V2[-2,j,k])/(2*dx)
                Wx[-1,j,k] = (W2[0,j,k] - W2[-2,j,k])/(2*dx)
            else:
                Ux[i,j,k] = (U2[i+1,j,k] - U2[i-1,j,k])/(2*dx)
                Vx[i,j,k] = (V2[i+1,j,k] - V2[i-1,j,k])/(2*dx)
                Wx[i,j,k] = (W2[i+1,j,k] - W2[i-1,j,k])/(2*dx)    
                
            if j==0:
                Uy[i,0,k] = (U2[i,1,k] - U2[i,-1,k])/(2*dx)
                Vy[i,0,k] = (V2[i,1,k] - V2[i,-1,k])/(2*dx)
                Wy[i,0,k] = (W2[i,1,k] - W2[i,-1,k])/(2*dx)
            elif j==(N-2):
                Uy[i,-1,k] = (U2[i,0,k] - U2[i,-2,k])/(2*dx)
                Vy[i,-1,k] = (V2[i,0,k] - V2[i,-2,k])/(2*dx)
                Wy[i,-1,k] = (W2[i,0,k] - W2[i,-2,k])/(2*dx)
            else:
                Uy[i,j,k] = (U2[i,j+1,k] - U2[i,j-1,k])/(2*dx)
                Vy[i,j,k] = (V2[i,j+1,k] - V2[i,j-1,k])/(2*dx)
                Wy[i,j,k] = (W2[i,j+1,k] - W2[i,j-1,k])/(2*dx)   
                
            if k==0:
                Uz[i,j,0] = (U2[i,j,1] - U2[i,j,-1])/(2*dx)
                Vz[i,j,0] = (V2[i,j,1] - V2[i,j,-1])/(2*dx)
                Wz[i,j,0] = (W2[i,j,1] - W2[i,j,-1])/(2*dx)
            elif k==(N-2):
                Uz[i,j,-1] = (U2[i,j,0] - U2[i,j,-2])/(2*dx)
                Vz[i,j,-1] = (V2[i,j,0] - V2[i,j,-2])/(2*dx)
                Wz[i,j,-1] = (W2[i,j,0] - W2[i,j,-2])/(2*dx)
            else:
                Uz[i,j,k] = (U2[i,j,k+1] - U2[i,j,k-1])/(2*dx)
                Vz[i,j,k] = (V2[i,j,k+1] - V2[i,j,k-1])/(2*dx)
                Wz[i,j,k] = (W2[i,j,k+1] - W2[i,j,k-1])/(2*dx) 

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


p2 = ptilde1


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
U1a = U1[:,:,0:42]
V1a = V1[:,:,0:42]
W1a = W1[:,:,0:42]
P1a = p1[:,:,0:42]

#Chunk2
U2a = U1[:,:,42:84]
V2a = V1[:,:,42:84]
W2a = W1[:,:,42:84]
P2a = p1[:,:,42:84]

#Chunk3
U3a = U1[:,:,84:126]
V3a = V1[:,:,84:126]
W3a = W1[:,:,84:126]
P3a = p1[:,:,84:126]

#Chunk4
U4a = U2[:,:,0:42]
V4a = V2[:,:,0:42]
W4a = W2[:,:,0:42]
P4a = p2[:,:,0:42]

#Chunk5
U5a = U2[:,:,42:84]
V5a = V2[:,:,42:84]
W5a = W2[:,:,42:84]
P5a = p2[:,:,42:84]

#Chunk6
U6a = U2[:,:,84:126]
V6a = V2[:,:,84:126]
W6a = W2[:,:,84:126]
P6a = p2[:,:,84:126]

totalX = 512
currentX = 128*6
totalOverlap = currentX - totalX


#%%

Ufinal = np.empty((totalX,N-1,42),dtype=np.double)
Vfinal = np.empty((totalX,N-1,42),dtype=np.double)
Wfinal = np.empty((totalX,N-1,42),dtype=np.double)
Pfinal = np.empty((totalX,N-1,42),dtype=np.double)


#%%

Ufinal[42:86,:,:] = U1a[42:86,:,:]
Vfinal[42:86,:,:] = V1a[42:86,:,:]
Wfinal[42:86,:,:] = W1a[42:86,:,:]
Pfinal[42:86,:,:] = P1a[42:86,:,:]

Ufinal[128:172,:,:] = U2a[42:86,:,:]
Vfinal[128:172,:,:] = V2a[42:86,:,:]
Wfinal[128:172,:,:] = W2a[42:86,:,:]
Pfinal[128:172,:,:] = P2a[42:86,:,:]

Ufinal[214:258,:,:] = U3a[42:86,:,:]
Vfinal[214:258,:,:] = V3a[42:86,:,:]
Wfinal[214:258,:,:] = W3a[42:86,:,:]
Pfinal[214:258,:,:] = P3a[42:86,:,:]

Ufinal[300:344,:,:] = U4a[42:86,:,:]
Vfinal[300:344,:,:] = V4a[42:86,:,:]
Wfinal[300:344,:,:] = W4a[42:86,:,:]
Pfinal[300:344,:,:] = P4a[42:86,:,:]

Ufinal[386:430,:,:] = U5a[42:86,:,:]
Vfinal[386:430,:,:] = V5a[42:86,:,:]
Wfinal[386:430,:,:] = W5a[42:86,:,:]
Pfinal[386:430,:,:] = P5a[42:86,:,:]

Ufinal[472:512,:,:] = U6a[42:82,:,:]
Vfinal[472:512,:,:] = V6a[42:82,:,:]
Wfinal[472:512,:,:] = W6a[42:82,:,:]
Pfinal[472:512,:,:] = P6a[42:82,:,:]


for i in range(0,42):
    #beta = 1 - np.cos((np.pi/2.0)*float(i)/41.0)
    beta = float(i)/41.0
    theta = (np.pi/2.0)*beta
    Ufinal[i,:,:] = np.cos(theta)*U6a[82+i,:,:] + np.sin(theta)*U1a[i,:,:]
    Vfinal[i,:,:] = np.cos(theta)*V6a[82+i,:,:] + np.sin(theta)*V1a[i,:,:]
    Wfinal[i,:,:] = np.cos(theta)*W6a[82+i,:,:] + np.sin(theta)*W1a[i,:,:]
    Pfinal[i,:,:] = np.cos(theta)*P6a[82+i,:,:] + np.sin(theta)*P1a[i,:,:]
    
    Ufinal[86+i,:,:] = np.cos(theta)*U1a[86+i,:,:] + np.sin(theta)*U2a[i,:,:]
    Vfinal[86+i,:,:] = np.cos(theta)*V1a[86+i,:,:] + np.sin(theta)*V2a[i,:,:]
    Wfinal[86+i,:,:] = np.cos(theta)*W1a[86+i,:,:] + np.sin(theta)*W2a[i,:,:]
    Pfinal[86+i,:,:] = np.cos(theta)*P1a[86+i,:,:] + np.sin(theta)*P2a[i,:,:]    

    Ufinal[172+i,:,:] = np.cos(theta)*U2a[86+i,:,:] + np.sin(theta)*U3a[i,:,:]
    Vfinal[172+i,:,:] = np.cos(theta)*V2a[86+i,:,:] + np.sin(theta)*V3a[i,:,:]
    Wfinal[172+i,:,:] = np.cos(theta)*W2a[86+i,:,:] + np.sin(theta)*W3a[i,:,:]
    Pfinal[172+i,:,:] = np.cos(theta)*P2a[86+i,:,:] + np.sin(theta)*P3a[i,:,:]    

    Ufinal[258+i,:,:] = np.cos(theta)*U3a[86+i,:,:] + np.sin(theta)*U4a[i,:,:]
    Vfinal[258+i,:,:] = np.cos(theta)*V3a[86+i,:,:] + np.sin(theta)*V4a[i,:,:]
    Wfinal[258+i,:,:] = np.cos(theta)*W3a[86+i,:,:] + np.sin(theta)*W4a[i,:,:]
    Pfinal[258+i,:,:] = np.cos(theta)*P3a[86+i,:,:] + np.sin(theta)*P4a[i,:,:]
    
    Ufinal[344+i,:,:] = np.cos(theta)*U4a[86+i,:,:] + np.sin(theta)*U5a[i,:,:]
    Vfinal[344+i,:,:] = np.cos(theta)*V4a[86+i,:,:] + np.sin(theta)*V5a[i,:,:]
    Wfinal[344+i,:,:] = np.cos(theta)*W4a[86+i,:,:] + np.sin(theta)*W5a[i,:,:]
    Pfinal[344+i,:,:] = np.cos(theta)*P4a[86+i,:,:] + np.sin(theta)*P5a[i,:,:]

    Ufinal[430+i,:,:] = np.cos(theta)*U5a[86+i,:,:] + np.sin(theta)*U6a[i,:,:]
    Vfinal[430+i,:,:] = np.cos(theta)*V5a[86+i,:,:] + np.sin(theta)*V6a[i,:,:]
    Wfinal[430+i,:,:] = np.cos(theta)*W5a[86+i,:,:] + np.sin(theta)*W6a[i,:,:]
    Pfinal[430+i,:,:] = np.cos(theta)*P5a[86+i,:,:] + np.sin(theta)*P6a[i,:,:]  



#%%

f = open('U_uprime1_N128_k8_512x128x42.dat','w');
g = open('V_uprime1_N128_k8_512x128x42.dat','w');
h = open('W_uprime1_N128_k8_512x128x42.dat','w');
pp = open('P_uprime1_N128_k8_512x128x42.dat','w')
for k in range(0,42):
    for j in range(0,128):
        for i in range(0,512):
            f.write("".join([str(Ufinal[i,j,k]), "\n"]))
            g.write("".join([str(Vfinal[i,j,k]), "\n"]))
            h.write("".join([str(Wfinal[i,j,k]), "\n"]))
            pp.write("".join([str(Pfinal[i,j,k]), "\n"]))
f.close()
g.close()
h.close()
pp.close()