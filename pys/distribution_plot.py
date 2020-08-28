#!/usr/bin/python3
import numpy as np


# definitions
def read_data(indir): # read .plt data
	f= open(indir,'r')
	first_line= f.readline()
	first_line= f.readline()
	rad=[]
	tht=[]
	zet=[]
	eng=[]
	pch=[]
	first_line= f.readline()
	for line in f:
		l= line.split()
		rad.append(float(l[0]))
		tht.append(float(l[1]))
		zet.append(float(l[2]))
		eng.append(float(l[3]))
		pch.append(float(l[4]))
	f.close()
	return rad,tht,zet,eng,pch

def gaussian(x,t,dff,xm): # gaussian function
	return np.exp(-(x -xm)**2/4/dff/t)/np.sqrt(4*np.pi*dff*t)

def maxwellian(en): # maxwellian function
	kt=1
	return 2*np.sqrt(en/np.pi)*(1/kt)**1.5*np.exp(-en/kt)



def format_func(value, tick_number): # find number of multiples of pi/2
    N= int(np.round(2*value/np.pi))
    if N==0:
        return "0"
    elif N==1:
        return r"$\pi/2$"
    elif N==2:
        return r"$\pi$"
    elif N==-1:
        return r"$-\pi/2$"
    elif N==-2:
        return r"$-\pi$"
    elif N%2>0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

# argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-file',required=True,
	help='input file')
parser.add_argument('-o','--output-file',
	help='output file')
parser.add_argument('-t','--type',
	help='plot type \'pos\' or \'vel\'',required=True)
parser.add_argument('-s','--show-plot',action="store_true",
	help='show_plot')
args = parser.parse_args()


# plot
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import pandas as pd
import seaborn as sns
sns.set()
plt.style.use('seaborn-whitegrid')
params = {'axes.labelsize' : 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large'}
pylab.rcParams.update(params)
rad,tht,zet,eng,pch= read_data(args.input_file)

if args.type=='pos': # position space distributions
	print("distribution plot pos...")
	fig, ax= plt.subplots(3,1,figsize=(6,6))
	from scipy.stats import norm
	rad_plot= sns.distplot(rad,ax=ax[0],
		fit=norm,kde=False)
	ax[0].set_xlabel('$r/a$')
	ax[0].set_ylabel('$f(r/a)$')
	ax[0].xaxis.set_major_locator(plt.MaxNLocator(7))
	ax[0].yaxis.set_major_locator(plt.MaxNLocator(4))
	ax[0].set(xlim= (0.2,0.8),ylim=(0,40))

	tht_plot= ax[1].hist(tht,bins=20,density=True,
		alpha=0.5)
	ax[1].set_xlabel('$\\theta$')
	ax[1].set_ylabel('$f(\\theta)$')
	ax[1].xaxis.set_major_locator(plt.MaxNLocator(8))
	ax[1].yaxis.set_major_locator(plt.MaxNLocator(4))
	ax[1].xaxis.set_major_locator(plt.MultipleLocator(np.pi/2))
	ax[1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
	ax[1].set(xlim=(-np.pi,np.pi),ylim=(0.0,0.25))

	zet_plot= ax[2].hist(zet,bins=20,density=True,
		alpha=0.5)
	ax[2].set_xlabel('$\zeta$')
	ax[2].set_ylabel('$f(\zeta)$')
	ax[2].xaxis.set_major_locator(plt.MaxNLocator(8))
	ax[2].yaxis.set_major_locator(plt.MaxNLocator(4))
	ax[2].xaxis.set_major_locator(plt.MultipleLocator(np.pi/2))
	ax[2].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
	ax[2].set(xlim=(0.0,2*np.pi),ylim=(0.0,0.25))

if args.type=='vel': # velocity space distributions
	print("distribution plot vel...")
	data= {'x': pch, 'y': eng}
	data = pd.DataFrame(data, columns=['x', 'y'])
	with sns.axes_style('white'):
		plot= sns.jointplot('x','y', data, kind='hex',
			xlim=(-1.0,1.0),ylim=(0.0,12.0))
		plot.set_axis_labels("$\lambda$","$\epsilon$(keV)")

plt.tight_layout()


# output
if args.output_file:
	plt.savefig(args.output_file)
if args.show_plot:
	plt.show()
