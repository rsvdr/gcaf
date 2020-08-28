#!/usr/bin/python3
import numpy as np


# definitions
def slice_data(pltfile): #read surf file
    x=[]
    y=[]
    f=open(pltfile,'r')
    for line in f:
    	L= line.split()
    	x.append(float(L[0]))
    	y.append(float(L[1]))
    f.close()
    return x,y


def triad_data(input,size): # triads from slices
    triad= []
    for i in range(0,size):
        surf_file= "/slice_"+str(i)+".plt"
        x= []
        y= []
        z= []
        x,z= slice_data(input+surf_file)
        y= [i*2*np.pi/size]*len(x)
        for j in range(len(x)):
            triad.append([x[j],y[j],z[j]])
    return triad


def vectors(triad): # make vectors from triads
    X= []
    Y= []
    Z= []
    for i in range(len(triad)):
        X.append(triad[i][0])
        Y.append(triad[i][1])
        Z.append(triad[i][2])
    return X,Y,Z


# argument parsing
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-folder',nargs='+',help='input folder',required=True)
parser.add_argument('-o','--output-file',help='output file',required=True)
parser.add_argument('-N','--number-of-files',
    type=int,help='number of files',required=True)
args = parser.parse_args()


# collect data
triad= triad_data(args.input_folder[0],args.number_of_files)
X,Y,Z= vectors(triad)


# write data
f=open(args.output_file,'w')
for i in range(len(X)):
    f.write(str(X[i])+"  "+str(Y[i])+"  "+str(Z[i])+'\n')
f.close()
