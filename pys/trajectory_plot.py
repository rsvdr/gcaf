#!/usr/bin/python3

def read_log(file, kwargs):
    """ read data file from gcafpp """
    import configparser
    config = configparser.ConfigParser()
    config.read(file)
    dum= config.get('field','minor_radius').split()
    kwargs['rmin']= float(dum[0])
    dum= config.get('field','major_radius').split()
    kwargs['rmaj']= float(dum[0])
    dum= config.get('particle','pitch_angle').split()
    kwargs['ptch']= dum[1]

# log data dictionary
log_data= {
    'rmin': 11,
    'rmaj': 111,
    'ptch': 1.1
    }

# argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-file',nargs='+',
    help='input file',required=True)
parser.add_argument('-t','--type',
	help='plot type \'xvst\' or \'rvsz\'',required=True)
parser.add_argument('-l','--log-file',nargs='+',
	help='log file info')
parser.add_argument('-o','--output-file',
	help='output file')
parser.add_argument('-s','--show-plot',action="store_true",
	help='show_plot')
args = parser.parse_args()
print("trajectory plot...")


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab
plt.style.use('seaborn-whitegrid')
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

if args.type=='xvst': # plot x vs t
	# read trajectory data
	tim=[]
	rad=[]
	tht=[]
	zet=[]
	eng=[]
	pch=[]
	f= open(args.input_file[0],'r')
	first_line= f.readline()
	for line in f:
		L= line.split()
		tim.append(float(L[0]))
		rad.append(float(L[1]))
		tht.append(float(L[2]))
		zet.append(float(L[3]))
		eng.append(float(L[4]))
		pch.append(float(L[5]))
	f.close()
	tim_scaled= [a*1000 for a in tim]

	fig, ax= plt.subplots(5,sharex=True,figsize=(6,8))
	for i in range(5):
		#ax[i].yaxis.set_major_locator(plt.MaxNLocator(3))
		ax[i].set(xlim=(min(tim_scaled),max(tim_scaled)))
	ax[0].plot(tim_scaled,rad)
	ax[0].set(ylabel='$r/a$',ylim=(0.0,1.0))
	ax[1].plot(tim_scaled,tht)
	ax[1].set(ylabel='$\\theta$',ylim=(-np.pi,np.pi))
	ax[2].plot(tim_scaled,zet)
	ax[2].set(ylabel='$\zeta$',ylim=(0.0,2*np.pi))
	ax[3].plot(tim_scaled,eng)
	ax[3].set(ylabel='$E(keV)$')
	ax[4].plot(tim_scaled,pch)
	ax[4].xaxis.set_major_locator(plt.MaxNLocator(5))
	ax[4].set(xlabel='$t(ms)$',ylabel='$\lambda$',ylim=(-1.1,1.1))
	plt.tight_layout()
	fig.subplots_adjust(hspace=0.1)

elif args.type=='rvsz': # plot R vs Z
    for k in range(len(args.input_file)):
        # read trajectory data
        x=[]
        y=[]
        z=[]
        f= open(args.input_file[k],'r')
        first_line= f.readline()
        for line in f:
            L= line.split()
            x.append(float(L[6]))
            y.append(float(L[7]))
            z.append(float(L[8]))
        f.close()
        R_cyl=[]
        for i in range(len(x)):
            R_cyl.append(np.sqrt(x[i]**2 +y[i]**2))
        Z_cyl= [a for a in z]

        if args.log_file: # if cfg file plot wall and mid rad
            read_log(args.log_file,log_data)
            fig, ax= plt.subplots(figsize=(6,6))
            tht_dum= np.arange(0.0,2*np.pi,.01)
            rw= [log_data['rmaj'] \
                +log_data['rmin']*np.sin(a) for a in tht_dum]
            zw= [log_data['rmin']*np.cos(a) for a in tht_dum]
            rm= [log_data['rmaj'] \
                +0.5*log_data['rmin']*np.sin(a) for a in tht_dum]
            zm= [0.5*log_data['rmin']*np.cos(a) for a in tht_dum]
            ax.plot(R_cyl,Z_cyl,label="$\lambda=$ "+log_data['ptch'])
            ax.plot(rw,zw,c='k',lw=2)
            ax.plot(rm,zm,c='k',ls='--',lw=2)
            plt.legend(framealpha=1,frameon=True)
        else:
            fig, ax= plt.subplots(figsize=(6,6))
        ax.plot(R_cyl,Z_cyl)
        ax.set_xlabel('$R(cm)$')
        ax.set_ylabel('$Z(cm)$')
        plt.tight_layout()

else:
    print("Incorrect plot type: "+args.type)

# output
if args.output_file:
    plt.savefig(args.output_file)
if args.show_plot:
    plt.show()
