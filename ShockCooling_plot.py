#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 07:36:49 2022

Plot the data generated using SchoCooling_full.py

@author: ben
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import interp1d

#Get the loss function 
filename='ShockCooling_full_data.h5'
f=h5py.File(filename, "r")
shocks={}
for var in f.keys():
	print(var)
	shocks[var]=np.array(f[var])

Va2=1.0
Cs2=5.0/3.0*shocks['beta']*Va2/2.0
Vs2=0.5*(Cs2+Va2-np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])*np.cos(shocks['theta'])))
Vf2=0.5*(Cs2+Va2+np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])*np.cos(shocks['theta'])))

print(Cs2,Va2,np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])**2))

print('beta=',shocks['beta'],'theta=',shocks['theta'])

print(np.sqrt(Cs2),np.sqrt(Vs2),np.sqrt(Va2),np.sqrt(Vf2))

#Get the Hau-Sonnerup shock solutions
#a2dhau=np.linspace(0,2,100)
#a2uhau=

adcool=np.where(shocks['a2d'])

fig, ax = plt.subplots(figsize=(9.7, 6),dpi=300)  
ax.plot(shocks['a2d'],shocks['a2u'],'k.', markersize=4)
ax.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax.set_ylim(0,4)
ax.set_xlim(0,2)
ax.set_ylabel('$A^{u2}$',fontsize=14)
ax.set_xlabel('$A^{d2}$',fontsize=14)
plt.show()
#plt.savefig('ShockCooling_plot.png',dpi=300,bbox_inches='tight')
