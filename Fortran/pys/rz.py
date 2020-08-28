#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import funs
import sys


arg= sys.argv
com= funs.dirs(arg)
indir= com[0]
outdir= com[1]
pltname= com[2]
logfile= com[3]
number= com[6]
if len(logfile)>0:
	vec= funs.readlog(logfile)
	rmin= vec[0]
	rmaj= vec[1]
	ptch= vec[2]
	
fig, ax= plt.subplots(figsize=(6,6))
if number==1:
	if len(logfile)>0:
		funs.xyline(ax,indir,8,9,ptch,'','r')
	else:
		funs.xyline(ax,indir,8,9,'','','r')
if number>1:
	for i in range(0,number):
		#vec= funs.readlog(logfile[:-4]+'_'+str(i)+'.log')
		#ptch= vec[2]
		funs.xyline(ax,indir[:-4]+'_'+str(i)+'.plt',8,9,'','','r')
		print indir[:-4]+'_'+str(i)+'.plt'
if len(logfile)>0:
	tht= np.arange(0.0,2*np.pi,.01)
	rw= [rmaj+rmin*np.sin(a) for a in tht]
	zw= [rmin*np.cos(a) for a in tht]
	rm= [rmaj+0.5*rmin*np.sin(a) for a in tht]
	zm= [0.5*rmin*np.cos(a) for a in tht]
	#rp= [rmaj+0.6*rmin*np.sin(a) for a in tht]
	#zp= [0.6*rmin*np.cos(a) for a in tht]
	ax.plot(rw,zw,color='k',lw=2)
	ax.plot(rm,zm,color='k',ls='--',lw=2)
	#ax.plot(rp,zp,color='b',ls='--',lw=2)
	ax.set_xlim(rmaj-rmin,rmaj+rmin)
	ax.set_ylim(-rmin,rmin)
	plt.legend(title='pitch angle')

funs.mylabels(ax,'$R$(m)','$Z$(m)')
funs.myformat(ax)
funs.savemyfig(outdir,pltname)
plt.show()
