# renorm method
import numpy as np
import pandas as pd

R = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
E=2.5
H=( np.dot(R,R)-2*E*R+E*E*np.eye(4) )

O = np.eye(4)
dt = 0.001
N = 10000
rcdxp=np.zeros((4,N))
ang=0.7

x=np.array([1,2,3,4]/np.sqrt(30))*np.cos(ang)
p=np.array([1,2,3,4]/np.sqrt(30))*np.sin(ang)
for i in range(N):
    x += - np.dot(H,x)*dt
    p += - np.dot(H,p)*dt
    scl = np.sqrt(sum(x**2)+sum(p**2))
    x = x/scl
    p = p/scl
    rcdxp[:,i]=x**2+p**2
aaa=pd.DataFrame(rcdxp.T)
aaa.to_csv('arbi_E%.1f.csv'%E,header=None,sep=' ')
