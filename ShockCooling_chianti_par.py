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
from scipy.interpolate import PchipInterpolator
import multiprocessing

#Get the loss function 
#Chianti loss function
filename='lossfunc_photo_scott.h5'
with h5py.File(filename, "r") as f: 
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    tlf = list(f['temperature'])
    lf=list(f['rad_loss'])/np.max(f['rad_loss'])

nropoints=1000000

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

def radJumpSol(nr):
	rmax=500.0
	RhoJump=nr*(rmax-1)/(nropoints+1)+1
	print(RhoJump)
	a2d=[]
	a2u=[]
	errarr=[]
	dT=[]
	for i in range(1,nelements):
		#Hau-Sonnerup solution
		a2dhau=float(i)*stepsize
		a2uhau=RhoJump*a2dhau #Why doesnt this work?
		if a2uhau <=4:
			#Find the potential solutions to the loss jump
			solarr=lossSol(a2uhau/a2dhau, lffunc(T0), lfint)
			for j in solarr:
				LossTjump=tlfint[j]/T0
				AlgTjump=tempjump(a2dhau,a2uhau,beta,theta)
		
		        #find the temperatre residule
				if np.abs(LossTjump-AlgTjump) <= 1.0e-2:
					a2d.append(a2dhau);a2u.append(a2uhau);errarr.append(LossTjump-AlgTjump);dT.append(LossTjump)
					
	return(a2d,a2u,errarr,dT)

beta=0.1
theta=3.14/8.0
#T0=tlf[np.argmax(lf)]
T0=1.0e6
gamma=5.0/3.0

nelements=10001
admax=2.0

a2d=[]
a2u=[]#np.empty(nelements)
errarr=[]
dT=[]

a2uhau=np.empty(nelements)
a2dhau=np.empty(nelements)

stepsize=admax/nelements

#Get a higher res loss function
tlfint=np.logspace(1,9,10000)
#lffunc=interp1d(tlf,lf,kind='cubic')
lffunc=PchipInterpolator(tlf,lf)
lfint=lffunc(tlfint)

#rarr=np.linspace(1.01,8.01,1001)
#rarr=np.linspace(1.01,10.01,101)
#rmax=8.01

pool = multiprocessing.Pool(30)
a2d,a2u,errarr,dT=zip(*pool.map(radJumpSol, range(0,nropoints)))
pool.close()

a2darr=[]
a2uarr=[]
errf=[]
dTarr=[]
for i in range(0,nropoints):
	for j in range(0,np.size(a2d[i])):
		a2darr.append(a2d[i][j])
		a2uarr.append(a2u[i][j])
		errf.append(errarr[i][j])
		dTarr.append(dT[i][j])
#stop
#for RhoJump in rarr:



#Save the data
hf = h5py.File('ShockCooling_chianti_data_T_1000000_par.h5', 'w')
hf.create_dataset('a2d', data=a2darr)
hf.create_dataset('a2u', data=a2uarr)
hf.create_dataset('errarr', data=errf)
hf.create_dataset('dT', data=dTarr)
hf.create_dataset('theta', data=theta)
hf.create_dataset('beta', data=beta)
#hf.create_dataset('rmax',data=np.max(rarr))
#hf.create_dataset('rmin',data=np.min(rarr))
hf.close()

#fig, ax = plt.subplots(figsize=(9.7, 6),dpi=300)  
#ax.plot(a2darr,a2uarr,'k.', markersize=4)
#ax.plot([0,2],[0,2])
#ax.set_ylim(0,4)
#ax.set_xlim(0,2)
#ax.set_ylabel('$A^{u2}$',fontsize=14)
#ax.set_xlabel('$A^{d2}$',fontsize=14)
#plt.savefig('ShockCooling_chianti_T_230000_par.png',dpi=300,bbox_inches='tight')
