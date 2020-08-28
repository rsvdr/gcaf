#!/usr/bin/python3
import numpy as np
import random

# definitions
def xyz_data(pltfile):
    f=open(pltfile,'r')
    x=[]
    y=[]
    z=[]
    first_line= f.readline()
    first_line= f.readline()
    first_line= f.readline()
    for line in f:
        L= line.split()
        x.append(float(L[0]))
        y.append(float(L[1]))
        z.append(float(L[2]))
    f.close()
    return x,y,z

def xyz_traj(pltfile):
    f=open(pltfile,'r')
    x=[]
    y=[]
    z=[]
    first_line= f.readline()
    first_line= f.readline()
    first_line= f.readline()
    for line in f:
        L= line.split()
        x.append(float(L[6]))
        y.append(float(L[7]))
        z.append(float(L[8]))
    f.close()
    return x,y,z


# argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-is','--input-surf',nargs='+',
    help='input surface data')
parser.add_argument('-id','--input-dist',nargs='+',
    help='input distribution data')
parser.add_argument('-it','--input-traj',nargs='+',
    help='input trajectory data')
parser.add_argument('-o','--output-file',
    help='output file')
parser.add_argument('-s','--show-plot',action="store_true",
    help='show_plot')
parser.add_argument('-ind','--index',type=int,
    help='plot index')
args = parser.parse_args()


# plot
from mayavi import mlab
if args.show_plot==False: # turn off screen
    mlab.options.offscreen = True
myfig= mlab.figure(1,fgcolor=(0,0,0),bgcolor=(1,1,1),size=(800,800))
if args.input_surf: #surface plot
    for i in range(len(args.input_surf)):
        X,Y,Z= xyz_traj(args.input_surf[i])
        mlab.plot3d(X,Y,Z,
            tube_radius=0.4,
            opacity=0.1,
            color=(0,0,1))

c_list=((random.uniform(0.0,1.0),
    random.uniform(0.0,1.0),
    random.uniform(0.0,1.0)))
if args.input_traj: #trajectory plot
    for i in range(len(args.input_traj)):
        X,Y,Z= xyz_traj(args.input_traj[i])
        mlab.points3d(X[args.index],Y[args.index],Z[args.index],
            #tube_radius=0.8,
            scale_factor=4,
            color=(1,0,0))

if args.input_dist: #distribution plot
    for i in range(len(args.input_dist)):
        X,Y,Z= xyz_data(args.input_dist[i])
        pts= mlab.points3d(X,Y,Z,
            scale_factor=4,
            color=(0,0.9,0))

# tokamak
rmaj= 150
rmin= 50
# z axis
t= np.linspace(-3*rmin,3*rmin,10)
x= [0 for a in t]
y= [0 for a in t]
z= [a for a in t]
mlab.plot3d(x,y,z,
    tube_radius=1,
    color=(0,0,0))
# magnetic axis
rmaj= 150
t= np.linspace(0.0,1.0,1000)
x= [rmaj*np.cos(2*np.pi*a) for a in t]
y= [rmaj*np.sin(2*np.pi*a) for a in t]
z= [0 for a in t]
mlab.plot3d(x,y,z,
    tube_radius=1,
    color=(0,0,0))
# wall
frq= 400
cut= -0.75
t= np.linspace(0.0,1.0,2000)
x= [(rmaj +rmin*np.cos(frq*np.pi*a))*np.cos(cut*2*np.pi*a) for a in t]
y= [(rmaj +rmin*np.cos(frq*np.pi*a))*np.sin(cut*2*np.pi*a) for a in t]
z= [rmin*np.sin(frq*np.pi*a) for a in t]
pts= mlab.plot3d(x,y,z,
    tube_radius=1,
    color=(0,0,0),
    opacity=0.05)

mlab.view(35,60,800)
#mlab.move(0,0,-1.5)

# output
if args.output_file:
    mlab.savefig(args.output_file)
if args.show_plot:
    mlab.show()
