#!/usr/bin/python

from scipy import stats
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np


def disp(ax,indir,n1,n2):
	f=open(indir,'r')
	first_line= f.readline()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
	f.close()
	slope,intercept,r_value,p_value,std_err= stats.linregress(x,y)
	fyd= [intercept+ slope*a for a in x]
	ax.plot(x,y,color='r',lw=2)
	if n2<4:
		ax.plot(x,fyd,'k',alpha=1.0,ls='--',label='$D=$'+str('%.2e'%slope))
	
	
def dirs(arg):
	indir= ''
	outdir= ''
	pltname= ''
	logfile= ''
	indir2= ''
	number= 1
	show= 'off'
	for i in range(len(arg)):
		if arg[i]== '-i':
			indir= arg[i+1]
		if arg[i]== '-o':
			outdir= arg[i+1]
		if arg[i]== '-n':
			pltname= arg[i+1]
		if arg[i]== '-l':
			logfile= arg[i+1]
		if arg[i]== '-i2':
			indir2= arg[i+1]
		if arg[i]== '-nmb':
			number= int(arg[i+1])
		if arg[i]== '-s':
			show= arg[i+1]
	com= [indir,outdir,pltname,logfile,indir2,show,number]
	return com


def readlog(logfile):
	rmin= 0.0
	rmaj= 0.0
	ptch= 0.0
	rad= 0.0
	vec= []
	f= open(logfile,'r')
	for line in f:
		L= line.split()
		if len(L)>0 and L[0]== 'minor_radius':
			rmin= float(L[1])
		if len(L)>0 and L[0]== 'major_radius':
			rmaj= float(L[1])
		if len(L)>0 and L[0]== 'pitch_angle':
			ptch= float(L[1])
		if len(L)>0 and L[0]== 'radius':
			rad= float(L[1])
		if len(L)>0 and L[0]== 'density':
			den= float(L[1])
	f.close()
	vec= [rmin,rmaj,ptch,rad,den]
	return vec


def savemyfig(outdir,pltname):
	if len(outdir)>0 and len(pltname)>0:
		plt.savefig(outdir+pltname+'.png')
	if len(outdir)>0 and len(pltname)==0:
		plt.savefig(outdir+'/traj.png')


def myformat(ax):
	ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
	ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	ax.xaxis.set_ticks_position('both')
	ax.yaxis.set_ticks_position('both')
	ax.tick_params(labelsize=14,direction='in',length=6,width=2)
	for axis in ['top','bottom','left','right']:
	  ax.spines[axis].set_linewidth(2)
	plt.tight_layout()


def mylabels(ax,xlab,ylab):
	ax.set_xlabel(xlab,labelpad=0,size=18)
	ax.set_ylabel(ylab,labelpad=0,size=18)


def xyerrbar(ax,pltfile,n1,n2,n3):
	f=open(pltfile,'r')
	first_line= f.readline()
	x=[]
	y=[]
	yer=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
		yer.append(6*float(L[n3]))
	f.close()
	ax.errorbar(x,y,yerr=yer,fmt='o',color='r',capsize=2)


def xyscat(ax,pltfile,n1,n2,mrk,clr,siz,al):
	f=open(pltfile,'r')
	first_line= f.readline()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
	f.close()
	ax.scatter(x,y,color=clr,alpha=al,marker=mrk,lw=0,s=siz)
	

def xyline(ax,pltfile,n1,n2,lab,mrk,clr):
	f=open(pltfile,'r')
	first_line= f.readline()
	l= first_line.split()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
	f.close()
	ax.set_xlim(min(x),max(x))
	ax.plot(x,y,label=lab,marker=mrk,ls='-',lw=2,color=clr)
	ax.plot(x[0],y[0],marker='o',c='r',ms=12)
	ax.plot(x[-1],y[-1],marker='X',c='r',ms=12)
	#ax.set_ylabel(str(l[n2+2]),size=14)

