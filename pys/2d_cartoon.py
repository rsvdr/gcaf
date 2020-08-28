#!/usr/bin/python3
import numpy as np


# definitions
def data(pltfile,n1,n2): # scatter data
	f=open(pltfile,'r')
	first_line= f.readline()
	first_line= f.readline()
	first_line= f.readline()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
	f.close()
	return x,y

def poin_dist(pltfile): # poincare distribution data
	f=open(pltfile,'r')
	first_line= f.readline()
	first_line= f.readline()
	first_line= f.readline()
	rad=[]
	tht=[]
	for line in f:
		L= line.split()
		if(float(L[2])<2*np.pi/100):
			rad.append(float(L[0]))
			tht.append(float(L[1]))
	f.close()
	return rad,tht



 # argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-id','--input-dist',
	help='input distribution')
parser.add_argument('-ip','--input-poin',
	help='input poincare')
parser.add_argument('-is','--input-flux-surf',
	help='input flux surface')
parser.add_argument('-o','--output-file',
	help='output file')
parser.add_argument('-s','--show-plot',action="store_true",
	help='show_plot')
args = parser.parse_args()
print("2d cartoon...")


# plot
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
plt.style.use('seaborn-whitegrid')
params = {'axes.labelsize' : 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)

fig, ax = plt.subplots(figsize=(12,6))
if args.input_dist:
	rad,tht= poin_dist(args.input_dist)
	#rad,tht= data(args.input_dist,0,1)
	ax.scatter(rad,tht,c='r',marker='o',s=12,alpha=1)
if args.input_poin:
	rad,tht= data(args.input_poin,1,2)
	ax.scatter(rad,tht,c='k',marker='.',s=5,alpha=0.5)
if args.input_flux_surf:
	rad,tht= data(args.input_flux_surf,1,0)
	ax.scatter(rad,tht,c='b',marker='x',s=12,alpha=1)
plt.ylim(-np.pi,np.pi)
#plt.xlim(0.0,0.6)
ax.set_xlabel("$r/a$")
ax.set_ylabel("$\\theta$")
plt.tight_layout()


# output
if args.output_file:
	plt.savefig(args.output_file)
if args.show_plot:
	plt.show()
