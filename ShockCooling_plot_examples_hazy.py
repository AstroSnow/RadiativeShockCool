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

#Define a figure float
fig, ((ax,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(9.7*1.2, 9.7*1.2),dpi=300)

#########################################################################################
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
 
line1=ax.plot(tlf,lf,'k',markersize=4,label='Ideal MHD')
ax.axvline(x=1000.0)
ax.axvline(x=10000.0)
ax.axvline(x=100000.0)
ax.axvline(x=1000000.0)
ax.axvline(x=10000000.0)
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('Normalised $\\Lambda (T)$',fontsize=14)
ax.set_xlabel('$T[K]$',fontsize=14)
#la=ax.text(0.05,3.7,'(a)')
#ax.legend()
#plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
#plt.show()

#########################################################################################
#Get the loss function 
filename='ShockCooling_hazy_data_T_1000_par.h5'
f=h5py.File(filename, "r")
shocks={}
for var in f.keys():
	#print(var)
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
 
line1=ax2.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax2.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax2.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax2.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax2.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax2.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax2.set_ylim(0,4)
ax2.set_xlim(0,2)
ax2.set_ylabel('$A^{u2}$',fontsize=14)
ax2.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
la=ax2.text(0.05,3.7,'(b)')

#########################################################################################
#Get the loss function 
filename='ShockCooling_hazy_data_T_10000.h5'
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
 
line1=ax3.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax3.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax3.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax3.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax3.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax3.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax3.set_ylim(0,4)
ax3.set_xlim(0,2)
ax3.set_ylabel('$A^{u2}$',fontsize=14)
ax3.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
la=ax3.text(0.05,3.7,'(c)')

#########################################################################################
#Get the loss function 
filename='ShockCooling_hazy_data_T_100000.h5'
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
 
line1=ax4.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax4.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax4.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax4.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax4.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax4.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax4.set_ylim(0,4)
ax4.set_xlim(0,2)
ax4.set_ylabel('$A^{u2}$',fontsize=14)
ax4.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
la=ax4.text(0.05,3.7,'(d)')

#########################################################################################
#Get the loss function 
filename='ShockCooling_hazy_data_T_1000000.h5'
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
 
line1=ax5.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax5.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax5.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax5.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax5.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax5.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax5.set_ylim(0,4)
ax5.set_xlim(0,2)
ax5.set_ylabel('$A^{u2}$',fontsize=14)
ax5.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
la=ax5.text(0.05,3.7,'(d)')

#########################################################################################
#Get the loss function 
filename='ShockCooling_hazy_data_T_10000000.h5'
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
 
line1=ax6.plot(a2dhau,a2uhau,'k.',markersize=4,label='Ideal MHD')
line2=ax6.plot(shocks['a2d'],shocks['a2u'],'g.', markersize=4,label='Isothermal')
line3=ax6.plot(shocks['a2d'][adcool[0]],shocks['a2u'][adcool[0]],'b.', markersize=4,label='Cooling')
line4=ax6.plot(shocks['a2d'][adheat[0]],shocks['a2u'][adheat[0]],'r.', markersize=4,label='Heating')

ax6.fill([0,2,2],[0,2,0], edgecolor='black', hatch='//',facecolor='w')
#ax.fill([0,0,np.max(a2uhau)/shocks['rmax']],[0,np.max(a2uhau),np.max(a2uhau)], edgecolor='black', hatch='//',facecolor='w')
#ax.fill_between(a2dhau, 1, where=a2dhau-a2uhau > 0, facecolor='green', alpha=.5)
ax6.plot([0,2],[0,2])
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
ax6.set_ylim(0,4)
ax6.set_xlim(0,2)
ax6.set_ylabel('$A^{u2}$',fontsize=14)
ax6.set_xlabel('$A^{d2}$',fontsize=14)
#ax.legend()
plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
la=ax6.text(0.05,3.7,'(d)')

#plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
plt.savefig('ShockCooling_plot_examles_hazy.png',dpi=300,bbox_inches='tight')
