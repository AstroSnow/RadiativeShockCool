#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 07:42:36 2023

@author: ben
"""
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import interp1d

#Get the loss function 
filename='lossfunc_photo_scott.h5'

with h5py.File(filename, "r") as f: 
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    tlf = list(f['temperature'])
    lf=list(f['rad_loss'])/np.max(f['rad_loss'])

loss2=np.zeros(np.size(lf))
loss3=np.zeros(np.size(lf))
for i in range(0,np.size(tlf)):
    loss2[i]=lf[i]
    loss3[i]=lf[i]
    if np.log10(tlf[i]) >= 5.5:
        loss2[i]=lf[i]*(1.0-np.tanh((np.log10(tlf[i])-5.5)/0.1)**2)
        loss3[i]=lf[i]*(1.0-np.tanh((np.log10(tlf[i])-5.5)/0.4)**2)
	
fig, ax = plt.subplots(figsize=(9.7, 6),dpi=300)
pmod=ax.plot(tlf,loss2,'orange',linewidth=3,label='Modified')
p=ax.plot(tlf,lf,linewidth=3,label='Cooling Curve')
p2=ax.plot([1.0e3,1.0e7],[0.25,0.25],'k--')
lMax=np.argmax(lf)
p3=ax.plot([tlf[lMax],tlf[lMax]],[0,1],'b--',label='Initial temperature')

#Find intersects
soall=np.zeros_like(lf)+0.25
s2=lf-soall
s3=s2[0:-2]*s2[1:-1]
pvar=0
for i in np.argwhere(s3 <=0):
	print(tlf[int(i)+1])
	if pvar==0:
		p4=ax.plot([tlf[int(i)+1],tlf[int(i)+1]],[0,1],'g--',label='r=2 solution temperature')
		pvar=1
	else:
		p5=ax.plot([tlf[int(i)+1],tlf[int(i)+1]],[0,1],'g--')
		
ax.legend(fontsize=14,loc='upper left')
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.set_xlim(1.0e3,1.0e7)
ax.set_ylabel('normalised $\Lambda(T)$',fontsize=14)
ax.set_xlabel('T [K]',fontsize=14)
plt.savefig('CoolingCurvePlot.png',dpi=300,bbox_inches='tight')