#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import setinfo

fflg=str(sys.argv[1]) # the filename to plot
clmn=int(sys.argv[2]) # the column to plot

myinfo = setinfo.Infolog('map.rc')

a=pd.read_csv(fflg,sep='\s+',header=None)
statname='null'
if not setinfo.is_number(a.values[0,0]):
    a=pd.read_csv(fflg,sep='\s+')
    statname=str(a.columns.values[clmn])

xs=a.values[:,0]*myinfo.dk['mapdtime']
if(clmn==0):
    ys=a.values[:,1:]
else:
    ys=a.values[:,clmn]

if(clmn==0):
    nclmn = len(ys[0,:])
    for i in range(nclmn):
        plt.plot(xs,ys[:,i],ls='--',c=(i/(nclmn-1),0,1-i/(nclmn-1)), label="%d-th state"%(i+1))
    plt.xlabel("time [a.u.]")
    plt.ylabel(r"population of %d states"%(nclmn))
    plt.legend(loc=1)
else:
    plt.plot(xs,ys,'r--')
    plt.xlabel("time [a.u.]")
    plt.ylabel(r"population on $\mid %d \rangle$ state"%clmn)
plt.title("Four spin-orbitals model")
plt.savefig(fflg+str(clmn)+'pop.png')

if(sys.argv[3]=='Y'):
	plt.show()

