#!/usr/bin/python3
import numpy as np
from scipy import stats

# definitions
def read_data(indir): # read .plt data
	f=open(indir,'r')
	first_line= f.readline()
	tim= []
	rad_dev= []
	rad_skw= []
	rad_fla= []
	for line in f:
		L= line.split()
		tim.append(float(L[0]))
		rad_dev.append(float(L[1]))
		rad_skw.append(float(L[4]))
		rad_fla.append(float(L[7]))
	f.close()
	return tim,rad_dev,rad_skw,rad_fla


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
	help='show_plot')
args = parser.parse_args()
print("moments plot...")

# plot
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib as mpl
plt.style.use('seaborn-whitegrid')
params = {'axes.labelsize' : 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large'}
pylab.rcParams.update(params)
fig, ax= plt.subplots(3,1,sharex=True,figsize=(6,8))

ls_list=['-','--',':']
label_list=[]
if args.legend_name:
	for i in range(len(args.legend_name)):
		label_list.append(args.legend_name[i])
else:
	for i in range(len(args.input_file)):
		label_list.append("fit")

for i in range(len(args.input_file)):
	tim,rad_dev,rad_skw,rad_fla= read_data(args.input_file[i])
	x_new= [a*1000 for a in tim]

	ax[0].plot(x_new,rad_dev, # data plot
		alpha=0.5)
		#c='b',alpha=0.5)
	ax[0].set_ylabel('$\langle \Delta r^2\\rangle$')

	polyfit= np.polyfit(x_new,rad_dev,2)
	fitting= np.poly1d(polyfit)
	x_fit= np.linspace(x_new[0],x_new[-1],1000)
	y_fit= fitting(x_fit)
	ax[0].plot(x_fit,y_fit, # fit plot
		c='k',ls=ls_list[i],lw=2,label=label_list[i])

	ax[1].plot(x_new,rad_skw, # data plot
		alpha=0.5)
		#c='b',alpha=0.5)

	polyfit= np.polyfit(x_new,rad_skw,8)
	fitting= np.poly1d(polyfit)
	x_fit= np.linspace(x_new[0],x_new[-1],1000)
	y_fit= fitting(x_fit)
	ax[1].plot(x_fit,y_fit, # fit plot
		c='k',ls=ls_list[i],lw=2,label=label_list[i])

	ax[2].plot(x_new,rad_fla, # data plot
		alpha=0.5)
		#c='b',alpha=0.5)

	polyfit= np.polyfit(x_new,rad_fla,8)
	fitting= np.poly1d(polyfit)
	x_fit= np.linspace(x_new[0],x_new[-1],1000)
	y_fit= fitting(x_fit)
	ax[2].plot(x_fit,y_fit, # fit plot
		c='k',ls=ls_list[i],lw=2,label=label_list[i])


#ax[0].yaxis.set_major_locator(plt.MaxNLocator(4))
ax[0].legend()

ax[1].axhline(0,c='k',lw=2,ls='-.',alpha=0.5,
	label='$\langle \Delta r^3/\sigma^3\\rangle= 0$')
ax[1].set(ylabel= '$\langle \Delta r^3/\sigma^3\\rangle$')
#ax[1].yaxis.set_major_locator(plt.MaxNLocator(4))
ax[1].legend()

ax[2].axhline(3,c='k',lw=2,ls='-.',alpha=0.5,
	label='$\langle \Delta r^4/\sigma^4\\rangle =3$')
ax[2].set(xlabel= "$t(ms)$",
	ylabel= '$\langle \Delta r^4/\sigma^4\\rangle$')
#ax[2].yaxis.set_major_locator(plt.MaxNLocator(4))
ax[2].xaxis.set_major_locator(plt.MaxNLocator(6))
ax[2].legend()

fig.tight_layout()
fig.subplots_adjust(hspace=0.1)


# output
if args.output_file:
	plt.savefig(args.output_file)
if args.show_plot:
	plt.show()
