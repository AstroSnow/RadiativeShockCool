#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 10:15:37 2023

@author: ben
"""
from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
import shockid
import h5py
from matplotlib import colors

#fname = '/media/ben/SnowM2/mhd_rad/'
fname = '/media/snow/SnowM2/mhd_rad/'

ShockJumpRead=True
ContextPlot=False

if ContextPlot==True:
	ds=pipread(fname,40) #Get the simulation data
	#Create a h5 file of the shock data
	shocks=shockid.restoreShocks(''.join((fname,'shocks.h5')))
	shocksf=shockid.shockFilter(shocks, 3)
	
	xg=ds['xgrid'][shocks['xs']:shocks['xe']]
	yg=ds['ygrid'][shocks['ys']:shocks['ye']]
	ro=ds['ro_p'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	pr=ds['pr_p'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	vx=ds['vx_p'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	vy=ds['vy_p'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	vz=ds['vz_p'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	bx=ds['bx'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	by=ds['by'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	bz=ds['bz'][shocks['ys']:shocks['ye'],shocks['xs']:shocks['xe']]
	
	ds2={}
	ds2['xgrid']=xg
	ds2['ygrid']=yg
	ds2['ro']=ro
	ds2['vx']=vx
	ds2['vy']=vy
	ds2['vz']=vz
	ds2['bx']=bx
	ds2['by']=by
	ds2['bz']=bz
	ds2['pr']=pr
	
	
	fig, ax = plt.subplots(figsize=(9, 6))
	plt.contourf(np.log10(ro),levels=101,cmap='Greys')
	sl=ax.plot(shocksf['slow'][:,1],shocksf['slow'][:,0],color='r',linestyle='',marker='.',markersize=2.8)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	plt.savefig('shockCoolingTest_context.png',dpi=300,bbox_inches='tight')


if ShockJumpRead == True:
	with h5py.File(''.join((fname,'ShockJumpData.h5')), "r") as f:
		# Print all root level object names (aka keys) 
		# these can be group or dataset names 
		print("Keys: %s" % f.keys())
		# get first object name/key; may or may NOT be a group
		a_group_key = list(f.keys())[0]

	# get the object type for a_group_key: usually group or dataset
		print(type(f[a_group_key])) 
	
	    # If a_group_key is a group name, 
	    # this gets the object names in the group and returns as a list
		Tjarr = np.array(f['Tjarr'])
		rojarr=np.array(f['rojarr'])
		prjarr=np.array(f['prjarr'])
		
		#Only have the compression solutions
		rojarr=rojarr[np.where(rojarr>1.0)]
		Tjarr=Tjarr[np.where(rojarr>1.0)]
		prjarr=prjarr[np.where(rojarr>1.0)]
else:	
	Tjarr=[]
	T2jarr=[]
	rojarr=[]
	prjarr=[]
	for tpoint in range(0,np.size(shocksf['slow'][:,0])):
	
		print(tpoint,'/',np.size(shocksf['slow'][:,0]))
		shockLine=shockid.shockLine(shocksf['slow'][tpoint,:], ds2)
		ipre,ipos=shockid.prepostIndex(shockLine['normarr']['ro'],5)
		
		T=shockLine['normarr']['pr']/shockLine['normarr']['ro']
		
		Tj=T[ipos]/T[ipre]
		T2j=(T[ipos]-T[ipre])/T[ipre]
		roj=shockLine['normarr']['ro'][ipos]/shockLine['normarr']['ro'][ipre]
		prj=shockLine['normarr']['pr'][ipos]/shockLine['normarr']['pr'][ipre]
		
		Tjarr.append(Tj)
		T2jarr.append(T2j)
		rojarr.append(roj)
		prjarr.append(prj)
	
	#Save the Data	
	hf = h5py.File(''.join((fname,'ShockJumpData.h5')), 'w')
	hf.create_dataset('Tjarr', data=Tjarr)
	hf.create_dataset('T2jarr', data=T2jarr)
	hf.create_dataset('rojarr', data=rojarr)
	hf.create_dataset('prjarr', data=prjarr)
	hf.close()
	stop
	
fig, ax = plt.subplots(figsize=(9, 6))
#plt.hist(np.log10(Tjarr),bins=np.linspace(-2,2,150))
Tjar=np.array(Tjarr)
#ax.fill([-2,-2,np.log10(0.9),np.log10(0.9)],[0,210,210,0],'b')
#ax.fill([np.log10(0.9),np.log10(0.9),np.log10(1.1),np.log10(1.1)],[0,210,210,0],'g')
#ax.fill([2,2,np.log10(1.1),np.log10(1.1)],[0,210,210,0],'r')
N, bins, patches=ax.hist(np.log10(Tjar),bins=np.linspace(-2,2,150),color='k')
for i in range(0,73):
    patches[i].set_facecolor('b')
for i in range(73,76):    
    patches[i].set_facecolor('g')
for i in range(76, len(patches)):
    patches[i].set_facecolor('r')
ax.text(-1.9,200,'Cooling shocks',color='b')
ax.text(1.3,200,'Heating shocks',color='r')
#ax.hist(np.log10(Tjar[np.where(Tjar<0.9)]),bins=np.linspace(-2,2,150))
#ax.hist(np.log10(Tjar[np.where(Tjar>1.1)]),bins=np.linspace(-2,2,150),color='r')
#ax.plot(np.log10([1.1,1.1]),[0,600],'r')
#ax.plot(np.log10([0.9,0.9]),[0,600],'b')
ax.set_xlim(-2,2)
ax.set_ylim(0,210)
ax.set_xlabel('$\log_{10}(T^d/T^u)$')
ax.set_ylabel('Counts')
plt.savefig('shockCoolingTest.png',dpi=300,bbox_inches='tight')

fig, ax = plt.subplots(figsize=(9, 6))
sp=plt.scatter(np.log10(rojarr),np.log10(Tjarr),c=np.log10(prjarr),norm=colors.CenteredNorm(),cmap='PRGn')
plt.colorbar(sp,label='$\log_{10}(P^d/P^u)$')
plt.plot([-0.5,2.0],[0,0],'k')
ax.text(1.25,-1.8,'Cooling shocks',color='b')
ax.text(1.25,2.0,'Heating shocks',color='r')
ax.set_xlim(0.0,1.6)
ax.set_ylabel('$\log_{10}(T^d/T^u)$')
ax.set_xlabel('$\log_{10}(\\rho^d/\\rho^u)$')
plt.savefig('shockCoolingPressureScatter.png',dpi=300,bbox_inches='tight')