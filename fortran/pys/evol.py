#!/usr/bin/python

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

fig, ((ax1,ax6),(ax2,ax7),(ax3,ax8),(ax4,ax9),(ax5,ax10))= plt.subplots(5,2,figsize=(12,8))

sz= 14
ax1.set_ylabel('$\psi_p$',size=sz)
ax2.set_ylabel('$\\theta$',size=sz)
ax3.set_ylabel('$\zeta$',size=sz)
ax4.set_ylabel('$P_\parallel/mc$',size=sz)
ax5.set_ylabel('$P_\perp/mc$',size=sz)
ax6.set_ylabel('$\lambda$',size=sz)
ax7.set_ylabel('$E(keV)$',size=sz)
ax8.set_ylabel('$P_\zeta$',size=sz)
ax8.set_xlabel('$t(s)$',size=sz)
ax9.set_axis_off()
ax10.set_axis_off()
n=0
for axi in [ax1,ax2,ax3,ax4,ax5,ax6,ax7]:
	n=n+1
	funs.xyline(axi,indir,0,n,'','','r')
	funs.myformat(axi)

#funs.savemyfig(outdir,pltname)
plt.show()
