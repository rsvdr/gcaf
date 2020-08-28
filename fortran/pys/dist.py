#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import funs
import sys


def gaussian(x,t,dff,xm):
	return np.exp(-(x -xm)**2/4/dff/t)/np.sqrt(4*np.pi*dff*t)
	
	
def maxwell(en):
	kt=1
	return 2*np.sqrt(en/np.pi)*(1/kt)**1.5*np.exp(-en/kt)


def flat(x):
	return 0.5


arg= sys.argv
vec =funs.dirs(arg)
indir= vec[0]
outdir= vec[1]
pltname= vec[2]
logfile= vec[3]
show= vec[5]

if len(logfile)>0:
	vec =funs.readlog(logfile)
	rmin= vec[0]
	rad= vec[3]
	den= vec[4]

f= open(indir,'r')
first_line= f.readline()
first_line= f.readline()
l= first_line.split()
t= float(l[0])
dff= float(l[1])
x=[]
y=[]
z=[]
first_line= f.readline()
for line in f:
	l= line.split()
	x.append(float(l[0]))
	y.append(float(l[1]))
	z.append(float(l[2]))
f.close()	
	
	
fig, (ax1,ax2,ax3)= plt.subplots(3,figsize=(8,8))
if len(logfile)>0:
	x_gaussian= np.arange(30.0,70.0,.1)
	y_gaussian= [gaussian(a,t,dff,rad) for a in x_gaussian]
	x_maxwell= np.arange(0,8,.01)
	y_maxwell= [maxwell(a) for a in x_maxwell]
	x_flat= np.arange(-1.0,1.0,0.01)
	y_flat= [flat(a) for a in x_flat]
sz=8
fs= 18
htype= 'bar'
setnorm= True
nbn= 40

ax1.hist(x,bins=nbn,histtype=htype,normed=setnorm)
if len(logfile)>0:
	ax1.plot(x_gaussian,y_gaussian,color='r',lw=2)
	ax1.set_xlim(0.0,rmin)
ax1.set_xlabel('$r$',fontsize=fs)
ax1.set_ylabel('$f(r)$',fontsize=fs)
#ax1.set_xlim(40,60)

ax2.hist(y,bins=nbn,histtype=htype,normed=setnorm)
if len(logfile)>0:
	ax2.plot(x_maxwell,y_maxwell,color='r',lw=2)
ax2.set_xlabel('$\mathcal{E}(keV)$',fontsize=fs)
ax2.set_ylabel('$f(\mathcal{E})$',fontsize=fs)

ax3.hist(z,bins=nbn,histtype=htype,normed=setnorm)
if len(logfile)>0:
	ax3.plot(x_flat,y_flat,color='r',lw=2)
ax3.set_xlim(-1.0,1.0)
ax3.set_xlabel('$\cos(\\alpha)$',fontsize=fs)
ax3.set_ylabel('$f(\mathcal{\lambda})$',fontsize=fs)

for axi in [ax1,ax2,ax3]:
	funs.myformat(axi)


plt.tight_layout()
if len(outdir)>0:
	ax1.set_title('n='+str(den),size=18)
	ax1.set_title(pltname,size=18)
	if len(pltname)>0:
		plt.savefig(outdir+pltname+'.png')
	if len(pltname)==0:
		plt.savefig(outdir+'dist.png')
if show=='on':
	plt.show()
