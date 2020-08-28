#!/usr/bin/python3
import numpy as np
import random as rn

# definitions
def xy_data(pltfile): #plot function
	f=open(pltfile,'r')
	first_line= f.readline()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[0]))
		y.append(float(L[1]))
	f.close()
	return x,y


# argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-file',nargs='+',required=True,
	help='input file')
parser.add_argument('-o','--output-file',
	help='output file')
parser.add_argument('-N','--legend-name',nargs='+',
	help='legend names')
parser.add_argument('-s','--show-plot',action="store_true",
	help='show plot')
args = parser.parse_args()


# plot
print("flux plot...")
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('seaborn-whitegrid')
import matplotlib.pylab as pylab
params = {'axes.labelsize' : 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)
from scipy.interpolate import UnivariateSpline

fig, ax= plt.subplots(1,2,figsize=(12,6))
ls_list=['-','--',':']
label_list=[]
if args.legend_name:
	for i in range(len(args.legend_name)):
		label_list.append(args.legend_name[i])
else:
	for i in range(len(args.input_file)):
		label_list.append("fit")

for i in range(len(args.input_file)):
	x,y= xy_data(args.input_file[i])
	x= [a*1000 for a in x]
	y= [a*100 for a in y]
	ax[0].plot(x,y, # data plot
		lw=1,alpha=0.5)
		#lw=1,c='b',alpha=0.5)

	z = np.polyfit(x,y,2)
	f = np.poly1d(z)
	fp = f.deriv()
	x_new = np.linspace(x[0],x[-1],1000)
	y_new = f(x_new)
	yp_new = fp(x_new)
	ax[0].plot(x_new,y_new, # fit plot
		lw=2,ls=ls_list[i],
		c='k',label=label_list[i])
	ax[1].plot(x_new,yp_new, # fit derivative plot
		lw=2,ls=ls_list[i],
		c='k',label=label_list[i])

#ax[0].set_xlim(0.0,0.4)
#ax[0].set_ylim(0,40)
ax[0].set(xlabel="$t$(ms)", ylabel="No. of particles (%)")
ax[0].legend()
ax[1].set(xlabel="$t$(ms)", ylabel="Particle flux")
ax[1].legend()
plt.tight_layout()

# output
if args.output_file:
	plt.savefig(args.output_file)
if args.show_plot:
	plt.show()
