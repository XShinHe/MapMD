# scaling method
import numpy as np
import pandas as pd

H = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
H2 = np.dot(H,H)
O = np.eye(4)
dt = 0.001
N = 10000
rcdxp=np.zeros((4,N))
ang=0.5

x=np.array([1,2,3,4]/np.sqrt(30))*np.cos(ang)
p=np.array([1,2,3,4]/np.sqrt(30))*np.sin(ang)

for i in range(N):
    cc = np.dot(x,x)+np.dot(p,p)
    #print('b ',cc)
    
    mH=( np.dot(x, np.dot(H,x)) + np.dot(p, np.dot(H,p)) )
    mH2=( np.dot(x, np.dot(H2,x)) + np.dot(p, np.dot(H2,p)) )
    E=mH/cc #- 1/dt + np.sqrt(mH**2/cc**2-mH2/cc+1/dt**2)
    
    x += - np.dot(H,x)*dt + E*x*dt
    p += - np.dot(H,p)*dt + E*p*dt
    
    cc = np.dot(x,x)+np.dot(p,p)
    #print('e ',cc)
    
    rcdxp[:,i]=x**2+p**2
    
aaa=pd.DataFrame(rcdxp.T)
aaa.to_csv('scaling.csv',header=None,sep=' ')
