#!/usr/bin/python

import pylab
import matplotlib.pyplot as plt
import numpy as np
import funs
import sys


arg= sys.argv
vec =funs.dirs(arg)
indir= vec[0]
outdir= vec[1]
pltname= vec[2]
logfile= vec[3]

n= int(pltname)
fig, ax= plt.subplots(figsize=(6,4))
funs.xyline(ax,indir,1,n,indir[-8:],'o','r')
funs.xyline(ax,outdir,1,n,outdir[-8:],'o','b')
#funs.xyline(ax,pltname,1,2,pltname[-8:],'o','g')
#funs.xyerrbar(ax,indir,0,1,4)
#ax.set_xscale('log')
#ax.set_yscale('log')

#funs.mylabels(ax,'Number of particles','D')

#funs.mylabels(ax,'n','D')
plt.legend()
plt.tight_layout()
#funs.myformat(ax)
plt.show()
