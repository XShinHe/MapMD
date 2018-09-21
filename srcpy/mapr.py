import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

class vecc:
    def __init__(self, **args):
        # initial size
        if 'vsize' in args:
            self.vsize = args['vsize']
        else:
            print('warning, loss of matrix size')
            exit()
        
        # initial dim: complex, quaternion, octonion
        if 'cdim' in args:
            self.cdim = args['cdim']
        else:
            # default
            self.cdim = 2
        self.vc=np.zeros((self.vsize, self.cdim))
        
        # initial mode
        if 'imod' in args:
            self.imod = args['imod']
        else:
            # default
            self.imod = 0
    
        # running mode
        if 'rmod' in args:
            self.rmod = args['rmod']
        else:
            # default
            self.rmod = 0
        
        # population mode
        if 'pmod' in args:
            self.pmod = args['pmod']
        else:
            # default
            self.pmod = 0
            
    def getkeys(self, mydict):
        self.N = mydict['Nstep']
        self.dt = mydict['dt']
        self.argdict = mydict
    
    def getH(self,H):
        self.H = np.copy(H)
        
    def getHs(self,Hd,Hnd):
        self.Hd = np.copy(Hd)
        self.Hnd = np.copy(Hnd)
        
    def iniz(self):
        # for [w,x,y,z] : popu on w, y
        # equivalently for [x, y, px, py]: popu on x, px
        phangle = self.argdict['phangle']
        istate = self.argdict['istate']
        
        if self.cdim == 2 and self.imod ==0 :
            self.vc[istate,0] = np.cos(phangle)
            self.vc[istate,1] = np.sin(phangle)
        elif( self.cdim == 4 and self.imod == 0 ):
            self.vc[istate,0] = np.cos(phangle)
            self.vc[istate,1] = 0
            self.vc[istate,2] = np.sin(phangle)
            self.vc[istate,3] = 0
        elif( self.cdim == 4 and self.imod == 1 ):
            self.vc[istate,0] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,1] = np.sin(phangle)/np.sqrt(2)
            self.vc[istate,2] = -np.sin(phangle)/np.sqrt(2)
            self.vc[istate,3] = np.cos(phangle)/np.sqrt(2)
        elif( self.cdim == 4 and self.imod == 2 ):
            self.vc[istate,0] = np.cos(phangle)*np.cos(phangle)
            self.vc[istate,1] = np.sin(phangle)*np.cos(phangle)
            self.vc[istate,2] = np.sin(phangle)*np.cos(phangle)
            self.vc[istate,3] = np.sin(phangle)*np.sin(phangle)
        else:
            print('iniz error')
            exit()
            
    def popu(self):
        self.pop = np.zeros(self.vsize)
        
        if ( self.pmod ==0 ):
            self.pop = np.sum(self.vc**2,axis=1)
        elif (self.pmod == 1):
            self.pop = (self.vc[:,0]+self.vc[:,3])**2 + (self.vc[:,1]-self.vc[:,2])**2
        elif (self.pmod == 2 ):
            self.pop = self.vc[:,0]*self.vc[:,3]-self.vc[:,1]*self.vc[:,2]
        else:
            print('popu loss')
            pass
        return self.pop
    
    def runmd(self):
        self.iniz()
        rcd = np.zeros((self.N, self.vsize))
        
        for i in range(self.N):
            HVdt = np.dot(self.H, self.vc)*self.dt
            HV1dt = np.dot(self.Hd, self.vc)*self.dt
            HV2dt = np.dot(self.Hnd, self.vc)*self.dt
            
            if(self.rmod == 0): # 2-component equation
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt*0.5
                self.vc[:,1] += - np.dot(self.H, self.vc[:,0])*self.dt
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt*0.5
            elif(self.rmod == 1): # 4-component equation: Model I
                self.vc[:,0] += HVdt[:,1]
                self.vc[:,1] -= HVdt[:,0]
                self.vc[:,2] += HVdt[:,3]
                self.vc[:,3] -= HVdt[:,2]
            elif(self.rmod == 2): # 4-component equation: Model my version
                self.vc[:,0] -= HVdt[:,2]
                self.vc[:,1] -= HVdt[:,0]
                self.vc[:,2] += HVdt[:,3]
                self.vc[:,3] += HVdt[:,1]
            elif(self.rmod == 3): # 4-component equation: Model III
                self.vc[:,0] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] - HV2dt[:,1]
                self.vc[:,1] += 0.5*HV1dt[:,0]+0.5*HV1dt[:,3] + HV2dt[:,0]
                self.vc[:,2] += -0.5*HV1dt[:,0]-0.5*HV1dt[:,3] - HV2dt[:,3]
                self.vc[:,3] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] + HV2dt[:,2]
            elif(self.rmod == 4): # 4-component equation: Model IV
                self.vc[:,0] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] + HV2dt[:,2]
                self.vc[:,1] += 0.5*HV1dt[:,0]+0.5*HV1dt[:,3] + HV2dt[:,3]
                self.vc[:,2] += -0.5*HV1dt[:,0]-0.5*HV1dt[:,3] - HV2dt[:,0]
                self.vc[:,3] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] - HV2dt[:,1]
            elif(self.rmod == 5): # 4-component equation: Model V
                self.vc[:,0] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] + 0.5*HV2dt[:,2] - 0.5*HV2dt[:,1]
                self.vc[:,1] += 0.5*HV1dt[:,0]+0.5*HV1dt[:,3] + 0.5*HV2dt[:,3] + 0.5*HV2dt[:,0]
                self.vc[:,2] += -0.5*HV1dt[:,0]-0.5*HV1dt[:,3] - 0.5*HV2dt[:,0] - 0.5*HV2dt[:,3] 
                self.vc[:,3] += 0.5*HV1dt[:,2]-0.5*HV1dt[:,1] - 0.5*HV2dt[:,1] + 0.5*HV2dt[:,2]
            elif(self.rmod ==6):
                self.vc[:,0] += -HV1dt[:,1] + HV2dt[:,2]
                self.vc[:,1] += HV1dt[:,0] + HV2dt[:,3]
                self.vc[:,2] += -HV1dt[:,3] - HV2dt[:,0]
                self.vc[:,3] += HV1dt[:,2] - HV2dt[:,1]
            rcd[i,:] = self.popu()
            
        popdat = pd.DataFrame(rcd)
        popdat.to_csv('pop.dat',header=None,sep=' ')
            

idict={'dt':0.00001, 'Nstep':1000000, 'phangle':0.5 ,'istate':0}
lamb = 0.2

Hd = np.array([[10,0,0],[0,7,0],[0,0,2]])
Hnd = lamb*np.array([[0,1,1],[1,0,1],[1,1,0]])
H = Hd+Hnd
Nsize = len(H)
Ndim = 4

mapvc = vecc(vsize=Nsize, cdim=Ndim, imod=0, pmod=0, rmod=0)
mapvc.getkeys(idict)
mapvc.getH(H)
mapvc.getHs(Hd,Hnd)
mapvc.runmd()

os.system('python3 ../py/pop.py pop.dat 1 Y')
