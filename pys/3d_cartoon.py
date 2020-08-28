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
        x.append(float(L[1]))
        y.append(float(L[2]))
        z.append(float(L[3]))
    f.close()
    # pos= np.where(np.abs(np.diff(y)) > 3.0)[0]
    # for i in range(len(pos)):
    #     x[pos[i]]= np.nan
    #     y[pos[i]]= np.nan
    #     z[pos[i]]= np.nan
    # pos= np.where(np.abs(np.diff(z)) > 5.5)[0]
    # for i in range(len(pos)):
    #     x[pos[i]]= np.nan
    #     y[pos[i]]= np.nan
    #     z[pos[i]]= np.nan
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
parser.add_argument('-o','--output-file',help='output file')
parser.add_argument('-s','--show-plot',action="store_true",help='show_plot')
args = parser.parse_args()


# plot
from mayavi import mlab
if args.show_plot==False: # turn off screen
    mlab.options.offscreen = True
myfig= mlab.figure(1,fgcolor=(0,0,0),bgcolor=(1,1,1),size=(800,800))
my_scale= 50
if args.input_surf: #surface plot
    for i in range(len(args.input_surf)):
        X,Y,Z= xyz_data(args.input_surf[i])
        for j in range(len(Z)):
            Z[j]*= my_scale
        pts= mlab.points3d(X,Y,Z,Z)
        mesh= mlab.pipeline.delaunay2d(pts)
        pts.remove()
        surf = mlab.pipeline.surface(mesh,
            representation='surface',
            colormap='plasma',
            opacity=1)

if args.input_dist: #distribution plot
    for i in range(len(args.input_dist)):
        X,Y,Z= xyz_data(args.input_dist[i])
        for j in range(len(Z)):
            X[j]*= my_scale
        pts= mlab.points3d(Y,Z,X,
            scale_factor=0.05,
            color=(random.uniform(0.0,1.0),
                random.uniform(0.0,1.0),
                random.uniform(0.0,1.0)))

if args.input_traj: #trajectory plot
    for i in range(len(args.input_traj)):
        X,Y,Z= xyz_traj(args.input_traj[i])
        for j in range(len(Z)):
            X[j]*= my_scale
        mlab.points3d(Y,Z,X,
            scale_factor=0.05,
            color=(random.uniform(0.0,1.0),
                random.uniform(0.0,1.0),
                random.uniform(0.0,1.0)))

mlab.outline(extent=[-np.pi,np.pi,0.0,2*np.pi,0.45*my_scale,0.55*my_scale])
mlab.axes(ranges=[-np.pi,np.pi,0.0,2*np.pi,0.45,0.55])
mlab.xlabel("tht")
mlab.ylabel("zet")
mlab.zlabel("rad")
mlab.view(90,80,22,focalpoint=(0,0,0.5*my_scale))
mlab.move(0,-0.5,-1.5)


# output
if args.output_file:
    mlab.savefig(args.output_file)
if args.show_plot:
    mlab.show()
