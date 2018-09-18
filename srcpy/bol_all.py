import numpy as np
import pandas as pd
import math

bset=[0.1*i for i in range(1,200)]
bsize=len(bset)
rcd=np.zeros((4,bsize))

for i in range(bsize):
    beta = bset[i]
    
    H = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
    O1 = np.eye(4)
    O2 = np.eye(4)
    Ox = np.eye(4)
    P = np.eye(4)
    
    N = 10000
    dt1 = beta/N
    dt2 = np.sqrt(beta)/N
    
    for j in range(N):
        O1 += -np.dot(H,O1)*dt1
        O2 += -2*(j+0.5)*dt2*np.dot(H,O2)*dt2
        #
        P = np.dot(P,H)*beta/(j+1)
        Ox += (-1)**(j+1)*P
    rcd[0,i] = beta    
    rcd[1,i] = np.trace(O1)
    rcd[2,i] = np.trace(O2)
    rcd[3,i] = np.trace(Ox)
a=pd.DataFrame(rcd)
a.to_csv("bol_all.txt")
