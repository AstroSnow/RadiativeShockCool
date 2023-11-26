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
from matplotlib.legend_handler import HandlerLine2D

#Get the loss function 
#filename='ShockCooling_full_data_T_max.h5'
filename='ShockCooling_full_data_T_max_beta_0.5_theta_025pi.h5'
f=h5py.File(filename, "r")
shocks={}
for var in f.keys():
	print(var)
	shocks[var]=np.array(f[var])
f.close()

Va2=1.0
Cs2=5.0/3.0*shocks['beta']*Va2/2.0
Vs2=0.5*(Cs2+Va2-np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])*np.cos(shocks['theta'])))
Vf2=0.5*(Cs2+Va2+np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])*np.cos(shocks['theta'])))

print(Cs2,Va2,np.sqrt((Va2+Cs2)**2 -4.0*Va2*Cs2*np.cos(shocks['theta'])**2))

print('beta=',shocks['beta'],'theta=',shocks['theta'])

print(np.sqrt(Cs2),np.sqrt(Vs2),np.sqrt(Va2),np.sqrt(Vf2))

#Get the Hau-Sonnerup shock solutions
nelements=1001
stepsize=np.max(shocks['a2d'])/nelements
a2dhau=np.zeros(nelements)
a2uhau=np.zeros(nelements)
gamma=5.0/3.0
for i in range(1,nelements):
    #Hau-Sonnerup solution
    a2dhau[i]=float(i)*stepsize
    a2uhau[i]=(a2dhau[i]* ((gamma-1)/gamma *((gamma+1.0)/(gamma-1.0)-np.tan(shocks['theta'])**2)*(a2dhau[i]-1.0)**2 
                           +np.tan(shocks['theta'])**2*((gamma-1.0)/gamma *a2dhau[i]-1.0)*(a2dhau[i]-2.0)) -shocks['beta']/(np.cos(shocks['theta'])**2)*(a2dhau[i]-1.0)**2)
    a2uhau[i]=a2uhau[i]/((gamma-1.0)/gamma*(a2dhau[i]-1.0)**2/(np.cos(shocks['theta'])**2)-a2dhau[i]*np.tan(shocks['theta'])**2*((gamma-1)/gamma*a2dhau[i]-1.0))

adcool=np.where(shocks['dT'] <=0.9 )
adheat=np.where(shocks['dT'] >=1.1 )

def update_prop(handle, orig):
    handle.update_from(orig)
    handle.set_marker("")
    handle.set_linestyle('-')

fig, ax = plt.subplots(figsize=(9.7, 6),dpi=300)  
line1=ax.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax.set_ylim(0,4)
ax.set_xlim(0,2)
ax.set_ylabel('$A^{u2}$',fontsize=14)
ax.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
#plt.show()
#plt.savefig('ShockCooling_plot_T_max.png',dpi=300,bbox_inches='tight')
plt.savefig('ShockCooling_full_data_T_max_beta_0.5_theta_025pi.png',dpi=300,bbox_inches='tight')

