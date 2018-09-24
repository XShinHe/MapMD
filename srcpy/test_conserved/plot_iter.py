# plot for iteration

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

flg = str(sys.argv[1])
idx = int(sys.argv[2])

clr=[(1,0,0),(0,1,0),(0,0,1),(1,0,1),(1,1,0),(0,1,1),(1,1,1)]
lbl=['qVVseperated','qVVgrouped','qVVcrossing','cSimplectric']
for i in range(idx):
    n=pd.read_csv('%s_%d.dat'%(flg,i),header=None,sep='\s+')
    t=n.values[:,0]*0.001
    x=n.values[:,1]
    plt.plot(t,x,'--',c=clr[i],label=lbl[i])
plt.xlabel(r'$t\;\;[a.u.]$')
plt.ylabel('%s [a.u.]'%flg)
plt.legend(loc=1)
plt.savefig('%s_all.png'%flg)
plt.show()
