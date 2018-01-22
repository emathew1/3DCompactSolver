# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 18:52:40 2018

@author: mat.mathews
"""

#%%

import numpy as np
import matplotlib.pyplot as plt


filename = './uAvg.out'
filenamerho = './rhoAvg.out'
dataIn = np.loadtxt(filename,dtype=np.double)
dataInRho = np.loadtxt(filenamerho,dtype=np.double)
delU = 0.6

#%%

Nt = np.shape(dataIn)[0]
Ny = np.shape(dataIn)[1]
t  = np.empty(Nt,dtype=np.double)
uAvg   = np.empty((Nt,Ny-1),dtype=np.double)
rhoAvg = np.empty((Nt,Ny-1),dtype=np.double)

for i in range(0,Nt):
    t[i] = dataIn[i,0]
    uAvg[i] = dataIn[i,1:]
    rhoAvg[i] = dataInRho[i,1:]
    

#%%
momThick = np.empty(Nt,dtype=np.double)
for i in range(0,Nt):
    momThick[i] = np.trapz(rhoAvg[i]*(-uAvg[i,:]+delU/2)*(uAvg[i,:]+delU/2))/delU**2
    
#%%

A3 = [0, 0.9923076923076923, 60.3174603174603, 1.7769230769230777, 114.28571428571428, 2.815384615384616,142.8571428571429, 3.4153846153846157,171.4285714285714, 3.9923076923076923,200, 4.546153846153846,230.15873015873018, 5.0769230769230775,258.7301587301587, 5.515384615384616,287.3015873015873, 5.953846153846153,314.2857142857143, 6.438461538461539,342.8571428571429, 6.876923076923077,395.2380952380953, 7.661538461538461,453.96825396825386, 8.584615384615384, 511.1111111111111, 9.392307692307693]
AS3 = [-1.5873015873015959, 0.9923076923076923,80.9523809523809, 2.1461538461538474,160.3174603174603, 3.5999999999999996,198.41269841269838, 4.292307692307693,238.09523809523807, 4.961538461538462,271.4285714285714, 5.561538461538461,309.5238095238095, 6.253846153846154,347.61904761904765, 7.015384615384615,384.1269841269841, 7.707692307692308]

A3t = A3[0::2]
A3d = A3[1::2]
AS3t = AS3[0::2]
AS3d = AS3[1::2]
#%%

plt.plot(t*delU,momThick/momThick[ls
                                  0],A3t,A3d,'ok-',AS3t,AS3d,'sb-')

#%%
y = np.linspace(0,Ny-1,Ny-1)
yy = 0.3*np.tanh(-(y-(Ny-1)/2)/2)