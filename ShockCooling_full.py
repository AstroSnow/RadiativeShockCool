#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 07:36:49 2022

radiative shock jumps for the example of 

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

def tempjump(a2d,a2u,beta,theta):
    tjump=a2d/a2u*(1.0+(2.0/(beta*(1.0+np.square(np.tan(theta)))))
                   *(a2u-a2d+0.5*np.square(np.tan(theta))*(1.0-np.square((a2u-1.0)/(a2d-1.0)))))
    return tjump

def lossSol(r,tu,lfint):
	soall=np.zeros_like(lfint)+1.0/r**2
	s2=lfint-soall
	s3=s2[0:-2]*s2[1:-1]
	res=np.argwhere(s3 <=0)+1	
	return(res)

beta=0.1
theta=3.14/8.0
T0=tlf[np.argmax(lf)]
gamma=5.0/3.0

nelements=5001
admax=2.0

a2d=[]
a2u=[]#np.empty(nelements)
errarr=[]
dT=[]

a2uhau=np.empty(nelements)
a2dhau=np.empty(nelements)

stepsize=admax/nelements

#Get a higher res loss function
tlfint=np.logspace(3,7,10000)
lffunc=interp1d(tlf,lf,kind='cubic')
lfint=lffunc(tlfint)

rarr=np.linspace(1.01,10.01,10001)

for RhoJump in rarr:
	print(RhoJump)
	for i in range(1,nelements):
		#Hau-Sonnerup solution
		a2dhau[i]=float(i)*stepsize
		a2uhau[i]=RhoJump*a2dhau[i]
		if a2uhau[i] <=4:
			#Find the potential solutions to the loss jump
			solarr=lossSol(a2uhau[i]/a2dhau[i], T0, lfint)
			for j in solarr:
				LossTjump=tlfint[j]/T0
				AlgTjump=tempjump(a2dhau[i],a2uhau[i],beta,theta)
		
		        #find the temperatre residule
				if np.abs(LossTjump-AlgTjump) <= 1.0e-2:
					a2d.append(a2dhau[i]);a2u.append(a2uhau[i]);errarr.append(LossTjump-AlgTjump);dT.append(LossTjump)


#Save the data
hf = h5py.File('ShockCooling_full_data_2.h5', 'w')
hf.create_dataset('a2d', data=a2d)
hf.create_dataset('a2u', data=a2u)
hf.create_dataset('errarr', data=errarr)
hf.create_dataset('dT', data=dT)
hf.create_dataset('theta', data=theta)
hf.create_dataset('beta', data=beta)
hf.create_dataset('rmax',data=np.max(rarr))
hf.create_dataset('rmin',data=np.min(rarr))
hf.close()

fig, ax = plt.subplots(figsize=(9.7, 6),dpi=300)  
ax.plot(a2d,a2u,'k.', markersize=4)
ax.plot([0,2],[0,2])
ax.set_ylim(0,4)
ax.set_xlim(0,2)
ax.set_ylabel('$A^{u2}$',fontsize=14)
ax.set_xlabel('$A^{d2}$',fontsize=14)
plt.savefig('ShockCooling_full.png',dpi=300,bbox_inches='tight')
