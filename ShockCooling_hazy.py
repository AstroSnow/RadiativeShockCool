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
from scipy.interpolate import interp1d,PchipInterpolator

#Get the loss function 
#HAZY loss function
tf=np.linspace(1,9,81)
lf=[-28.233,-28.056,-27.871,-27.684,-27.494,-27.312,-27.136,-26.963,-26.783,
	-26.582,-26.362,-26.135,-25.922,-25.734,-25.572,-25.437,-25.324,-25.228,
	-25.146,-25.072,-25.003,-24.941,-24.876,-24.805,-24.739,-24.664,-24.603,
	-24.545,-24.447,-24.274,-23.898,-23.115,-22.028,-21.828,-21.932,-21.923,
	-21.833,-21.691,-21.510,-21.339,-21.281,-21.305,-21.288,-21.273,-21.276,
	-21.451,-21.666,-21.719,-21.779,-21.814,-21.788,-21.800,-21.885,-22.071,
	-22.253,-22.387,-22.469,-22.501,-22.508,-22.504,-22.508,-22.554,-22.623,
	-22.662,-22.669,-22.655,-22.629,-22.595,-22.557,-22.516,-22.473,-22.429,
	-22.381,-22.332,-22.279,-22.222,-22.162,-22.098,-22.028,-21.951,-21.866]
tlf=10**np.array(tf)
lf=10**np.array(lf)
lf=lf/np.max(lf)


def tempjump(a2d,a2u,beta,theta):
    tjump=a2d/a2u*(1.0+(2.0/(beta*(1.0+np.square(np.tan(theta)))))
                   *(a2u-a2d+0.5*np.square(np.tan(theta))*(1.0-np.square((a2u-1.0)/(a2d-1.0)))))
    return tjump

def lossSol(r,lftu,lfint):
	#soall=np.zeros_like(lfint)+1.0/r**2
	soall=np.zeros_like(lfint)+r**2
	#s2=lfint-soall
	s2=lftu/lfint-soall
	s3=s2[0:-2]*s2[1:-1]
	res=np.argwhere(s3 <=0)+1	
	return(res)

beta=0.1
theta=3.14/8.0
#T0=tlf[np.argmax(lf)]
T0=1.0e3
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
tlfint=np.logspace(2,7,10000)
#lffunc=interp1d(tlf,lf,kind='cubic')
lffunc=PchipInterpolator(tlf,lf)
lfint=lffunc(tlfint)

rarr=np.linspace(1.01,10.01,10001)
#rarr=np.linspace(1.01,10.01,101)

for RhoJump in rarr:
	print(RhoJump)
	a2uhau=np.zeros(nelements)
	a2dhau=np.zeros(nelements)
	for i in range(1,nelements):
		#Hau-Sonnerup solution
		a2dhau[i]=float(i)*stepsize
		a2uhau[i]=RhoJump*a2dhau[i] #Why doesnt this work?
		if a2uhau[i] <=4:
			#Find the potential solutions to the loss jump
			solarr=lossSol(a2uhau[i]/a2dhau[i], lffunc(T0), lfint)
			for j in solarr:
				LossTjump=tlfint[j]/T0
				AlgTjump=tempjump(a2dhau[i],a2uhau[i],beta,theta)
		
		        #find the temperatre residule
				if np.abs(LossTjump-AlgTjump) <= 1.0e-2:
					a2d.append(a2dhau[i]);a2u.append(a2uhau[i]);errarr.append(LossTjump-AlgTjump);dT.append(LossTjump)


#Save the data
hf = h5py.File('ShockCooling_hazy_data_T_1000.h5', 'w')
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
plt.savefig('ShockCooling_hazy_T_1000.png',dpi=300,bbox_inches='tight')
