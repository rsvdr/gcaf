#!/usr/bin/python3
import numpy as np
import argparse


# definitions
def data(pltfile,n1,n2): # read data
	f=open(pltfile,'r')
	first_line= f.readline()
	x=[]
	y=[]
	for line in f:
		L= line.split()
		x.append(float(L[n1]))
		y.append(float(L[n2]))
	f.close()
	return x,y

def takeFirst(elem):
    return elem[0]

def sort_poin(pltfile,n1,n2): # read and sort data
	Triple =[]
	# read data and create an unsorted list with [k,x,y] elements
	f=open(pltfile,'r')
	first_line= f.readline()
	for line in f:
		L= line.split()
		Triple.append([int(L[0]),float(L[n1]),float(L[n2])])
	f.close()
	# create a sorted list of [k,x,y] elements
	sortedTriple = sorted(Triple,key=takeFirst)
	# create X,Y lists for each k
	X=[]
	Y=[]
	tmp_X=[]
	tmp_Y=[]
	n=0
	for i in range(len(sortedTriple)):
		if sortedTriple[i][0]>n:
			X.append(tmp_X)
			Y.append(tmp_Y)
			tmp_X=[]
			tmp_Y=[]
			n+= 1
		tmp_X.append(sortedTriple[i][1])
		tmp_Y.append(sortedTriple[i][2])
	return X,Y

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
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-file',nargs='+',
	help='input file',required=True)
parser.add_argument('-t','--type',help='plot type \'rt\' or \'rt_color\' or surf or grid',required=True)
parser.add_argument('-o','--output-file',help='output file')
parser.add_argument('-s','--show-plot',action="store_true",help='show_plot')
args = parser.parse_args()
print("poincare plot...")


# plot
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'axes.labelsize':'x-large',
		'xtick.labelsize':'large',
        'ytick.labelsize':'large'}
pylab.rcParams.update(params)
if args.type=='rt': # rad vs tht
	fig, ax = plt.subplots(figsize=(6,6))
	for i in range(len(args.input_file)):
		x,y= data(args.input_file[i],1,2)
		ax.scatter(x,y,
			s=1,marker='.',alpha=0.5)
	ax.set_xlabel('$r/a$',size=18)
	ax.set_ylabel('$\\theta$',size=18)
	plt.ylim(-np.pi,np.pi)
	plt.tight_layout()
elif args.type=='rt_color': # rad vs tht colors
	fig, ax = plt.subplots(figsize=(6,6))
	x,y= sort_poin(args.input_file[0],1,2)
	for j in range(len(x)):
		ax.scatter(x[j],y[j],marker='.',s=4)
	# ax.set_xlabel('$r/a$',size=18)
	# ax.set_ylabel('$\\theta$',size=18)
	ax.set(xlim=(0.0,1.0),ylim=(-np.pi,np.pi))
elif args.type=='surf': # single surface
	fig, ax = plt.subplots(figsize=(6,6))
	data(args.input_file[0],1,0)
	ax.set_xlabel('$r/a$',size=18)
	ax.set_ylabel('$\\theta$',size=18)
	plt.ylim(-np.pi,np.pi)
	style(ax)
	plt.tight_layout()
elif args.type=='grid': # grid
	#plt.style.use('classic')
	fig, ax = plt.subplots(4,4,sharex=True,sharey=True,figsize=(8,8))
	n=0
	for i in range(4):
		for j in range(4):
			x,y= data(args.input_file[n],1,2)
			ax[i,j].xaxis.set_major_locator(plt.MaxNLocator(4))
			ax[i,j].yaxis.set_major_locator(plt.MaxNLocator(6))
			ax[i,j].yaxis.set_major_locator(plt.MultipleLocator(np.pi/2))
			ax[i,j].yaxis.set_major_formatter(plt.FuncFormatter(format_func))
			ax[i,j].scatter(x,y,
				c='k',s=0.1,marker='.')
			ax[i,j].set(xlim=(0.455,0.545),ylim=(-np.pi,np.pi))
			ax[i,0].set(ylabel="$\\theta$")
			ax[3,j].set(xlabel="$r/a$")
			n+= 1
	# text= ["(a)","(b)"]
	# n= 0
	# for i in range(2):
	# 	ax[i,0].text(-0.2,2.8,text[n],
	# 		color='k',fontsize=18)
	# 	n+=1
	# text= ["(c)","(d)"]
	# n= 0
	# for i in range(0,2):
	# 	for j in range(1,2):
	# 		ax[i,j].text(-0.14,2.8,text[n],
	# 			color='k',fontsize=18)
	# 		n+=1
	plt.tight_layout()
	fig.subplots_adjust(hspace=0.2, wspace=0.2)
else:
	print("Invalid type: "+args.type)

# output
if args.output_file:
	plt.savefig(args.output_file)
if args.show_plot:
	plt.show()
