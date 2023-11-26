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
	
pmod=ax.plot(tlf,loss2,'orange',linewidth=3,label='Modified')
p=ax.plot(tlf,lf,linewidth=3,label='Cooling Curve')
p2=ax.plot([1.0e3,1.0e7],[0.25,0.25],'k--')
lMax=np.argmax(lf)
#p3=ax.plot([tlf[lMax],tlf[lMax]],[0,1],'b--',label='Initial temperature')
p3=ax.plot([tlf[lMax],tlf[lMax]],[0,1],'b--',label='$T^u$')

#Find intersects
soall=np.zeros_like(lf)+0.25
s2=lf-soall
s3=s2[0:-2]*s2[1:-1]
pvar=0
for i in np.argwhere(s3 <=0):
	print(tlf[int(i)+1])
	if pvar==0:
		#p4=ax.plot([tlf[int(i)+1],tlf[int(i)+1]],[0,1],'g--',label='r=2 solution temperature')
		p4=ax.plot([tlf[int(i)+1],tlf[int(i)+1]],[0,1],'g--',label='r=2 solutions')
		pvar=1
	else:
		p5=ax.plot([tlf[int(i)+1],tlf[int(i)+1]],[0,1],'g--')
		
ax.legend(loc='upper left')
ax.set_xscale('log')
ax.set_ylim(0,1)
ax.set_xlim(1.0e3,1.0e7)
ax.set_ylabel('normalised $\Lambda(T)$',fontsize=14)
ax.set_xlabel('T [K]',fontsize=14)

#ax.axvline(x=1.8e4)
#ax.axvline(x=5e4)
#ax.axvline(x=2.3e5)
#ax.axvline(x=2.5e5)
#ax.axvline(x=5.0e5)
#ax.plot([0,2],np.sqrt([Va2,Va2]))
#ax.plot([0,2],[Cs2,Cs2])
#ax.plot([0,2],np.sqrt([Vs2,Vs2]))
#ax.plot([0,2],np.sqrt([Vf2,Vf2]))
#la=ax.text(0.05,3.7,'(a)')
#ax.legend()
#plt.legend(handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)})
#plt.show()

#########################################################################################
#Get the loss function 
filename='ShockCooling_full_data_T_18000.h5'
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
filename='ShockCooling_full_data_T_50000.h5'
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
filename='ShockCooling_full_data_T_max.h5'
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
filename='ShockCooling_full_data_T_250000.h5'
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
filename='ShockCooling_full_data_T_500000.h5'
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
plt.savefig('ShockCooling_plot_examles_chianti.png',dpi=300,bbox_inches='tight')
