#!/usr/bin/python
from scipy.interpolate import interp1d
import numpy as np
import argparse


# definitions
def takeFirst(elem):
    return elem[0]

def takeSecond(elem):
    return elem[1]

def PartTraj(pltfile,nprt): #returns trajectory points of nprt particle
	#read data points and create a list of [k,x,y] elements
	Triple =[]
	f=open(pltfile,'r')
	first_line= f.readline()
	for line in f:
		L= line.split()
		Triple.append([int(L[0]),float(L[1]),float(L[2])])
	f.close()
	#sort the list of [k,x,y] elements with ascending k
	sortedTriple = sorted(Triple,key=takeFirst)
	#create X,Y lists for each k
	X=[]
	Y=[]
	tmp_X=[]
	tmp_Y=[]
	n=0
	for i in range(len(sortedTriple)): #particle n=0
		X.append(tmp_X)
		Y.append(tmp_Y)
		tmp_X.append(sortedTriple[i][1])
		tmp_Y.append(sortedTriple[i][2])
	for i in range(len(sortedTriple)):
		if sortedTriple[i][0]>n:
			X.append(tmp_X)
			Y.append(tmp_Y)
			tmp_X=[]
			tmp_Y=[]
			n+= 1
		tmp_X.append(sortedTriple[i][1])
		tmp_Y.append(sortedTriple[i][2])
	return X[nprt],Y[nprt]


# argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input-file',help='input file',required=True)
parser.add_argument('-o','--output-file',help='output file',required=True)
parser.add_argument('-n','--particle-tag',type=int,help='particle tag',required=True)
parser.add_argument('-N','--file-name',help='file name')
args = parser.parse_args()

#create a list of [X,Y] elements and sort it
X,Y= PartTraj(args.input_file,args.particle_tag)
rt= []
for i in range(len(X)):
	rt.append([X[i],Y[i]])
rt_sorted= sorted(rt,key=takeSecond)

#obtain r, t arrays and interpolate
r=[]
t=[]
for i in range(len(rt_sorted)):
    r.append(rt_sorted[i][0])
    t.append(rt_sorted[i][1])
r_i1d= interp1d(t,r,kind='nearest',fill_value="extrapolate")
t_new= np.linspace(-3.14159265, 3.14159265, num=100)

#write the interpolated data
if args.file_name:
    f=open(args.output_file+args.file_name,'w')
else:
    f=open(args.output_file+"/slice.plt",'w')
for i in range(len(t_new)):
    f.write(str(t_new[i])+"  "+str(r_i1d(t_new)[i]))
    f.write('\n')
f.close()
