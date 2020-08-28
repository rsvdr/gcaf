#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import funs
import sys


arg= sys.argv
vec =funs.dirs(arg)
indir= vec[0]
outdir= vec[1]
pltname= vec[2]
logfile= vec[3]
indir2= vec[4]
show= vec[5]


ptsz= 4
fig, ax = plt.subplots(figsize=(6,6))
if logfile=='rt':
	funs.xyscat(ax,indir,1,2,'.','k',ptsz,1)
	if len(indir2)>0:
		funs.xyscat(ax,indir2,1,2,'.','r',2*ptsz,1)
	#plt.xlim(.48,.52)
	plt.ylim(-np.pi,np.pi)
if logfile=='rz':
	funs.xyscat(ax,indir,3,4,'.','k',ptsz,1)
	if len(indir2)>0:
		funs.xyscat(ax,indir2,3,4,'.','r',2*ptsz,1)

funs.mylabels(ax,'$r/a$','$\\theta$')
funs.myformat(ax)

plt.title(pltname,size=18)
plt.tight_layout()
if len(outdir)>0 and len(pltname)>0:
	plt.savefig(outdir+pltname+'.png')
if len(outdir)>0 and len(pltname)==0:
	plt.savefig(outdir+'traj.png')
if show=='on':
	plt.show()
