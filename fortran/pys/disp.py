#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import funs
import sys


arg= sys.argv
vec =funs.dirs(arg)
indir= vec[0]
outdir= vec[1]
pltname= vec[2]
logfile = vec[3]
show = vec[5]

fig, ((ax1,ax4,ax7),(ax2,ax5,ax8),(ax3,ax6,ax9))= plt.subplots(3,3,figsize=(12,6))
funs.disp(ax1,indir,0,1)
mylabel= '$\langle \Delta r^2\\rangle$'
ax1.set_ylabel(mylabel,size=14)

funs.disp(ax2,indir,0,2)
mylabel= '$\langle \Delta E^2\\rangle$'
ax2.set_ylabel(mylabel,size=14)

funs.disp(ax3,indir,0,3)
mylabel= '$\langle \Delta P_{\zeta}^2\\rangle$'
ax3.set_ylabel(mylabel,size=14)
ax3.set_xlabel('$t$(s)',size=14)

funs.disp(ax4,indir,0,4)
mylabel= '$\langle \Delta r^3/\sigma^3\\rangle$'
ax4.set_ylabel(mylabel,size=14)
ax4.set_ylim(-10.,10.)
ax4.axhline(0,color='k',lw=1)

funs.disp(ax5,indir,0,5)
mylabel= '$\langle \Delta E^3/\sigma^3\\rangle$'
ax5.set_ylabel(mylabel,size=14)

funs.disp(ax6,indir,0,6)
mylabel= '$\langle \Delta P_{\zeta}^3/\sigma^3\\rangle$'
ax6.set_ylabel(mylabel,size=14)
ax6.set_xlabel('$t$(s)',size=14)
ax6.axhline(0,color='k',lw=1)

funs.disp(ax7,indir,0,7)
mylabel= '$\langle \Delta r^4/\sigma^4\\rangle$'
ax7.set_ylabel(mylabel,size=14)
ax7.set_ylim(-10.,10.)
ax7.axhline(3,color='k',lw=1)

funs.disp(ax8,indir,0,8)
mylabel= '$\langle \Delta E^4/\sigma^4\\rangle$'
ax8.set_ylabel(mylabel,size=14)

funs.disp(ax9,indir,0,9)
mylabel= '$\langle \Delta P_{\zeta}^4/\sigma^4\\rangle$'
ax9.set_ylabel(mylabel,size=14)
ax9.set_xlabel('$t$(s)',size=14)
ax9.set_ylim(-10.,10.)
ax9.axhline(3,color='k',lw=1)

plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],visible=False)
for axi in [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]:
	funs.myformat(axi)
for axi in [ax1,ax2,ax3]:
	axi.legend()

ax1.set_title(pltname,size=18)
plt.tight_layout()
if len(outdir)>0 and len(pltname)>0:
	plt.savefig(outdir+pltname+'.png')
if len(outdir)>0 and len(pltname)==0:
	plt.savefig(outdir+'disp.png')
if show=='on':
	plt.show()
