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

fname = '/media/ben/SnowM2/mhd_rad/'
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
plt.savefig('shockCoolingTest_context.png',dpi=300)

Tjarr=[]
rojarr=[]
for tpoint in range(0,np.size(shocksf['slow'][:,0])):

	print(tpoint)
	shockLine=shockid.shockLine(shocksf['slow'][tpoint,:], ds2)
	ipre,ipos=shockid.prepostIndex(shockLine['normarr']['ro'],5)
	
	T=shockLine['normarr']['pr']/shockLine['normarr']['ro']
	
	Tj=T[ipos]/T[ipre]
	roj=shockLine['normarr']['ro'][ipos]/shockLine['normarr']['ro'][ipre]
	
	Tjarr.append(Tj)
	rojarr.append(roj)
	
fig, ax = plt.subplots(figsize=(9, 6))
ax.hist(Tjarr,bins=np.linspace(0,14,150))
ax.plot([1,1],[0,600],'r')
ax.set_xlim(0,15)
ax.set_ylim(0,600)
ax.set_xlabel('$T^d/T^u$')
ax.set_ylabel('Counts')
plt.savefig('shockCoolingTest.png',dpi=300)
