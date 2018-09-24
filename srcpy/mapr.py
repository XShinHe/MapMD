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
        
        if self.cdim == 2 and self.imod ==0 :  # 2-component normal initial
            self.vc[istate,0] = np.cos(phangle)
            self.vc[istate,1] = np.sin(phangle)
            self.pmod = 0
        elif( self.cdim == 4 and self.imod == 0 ): # 4-component normal initial # also other kind of initial methods
            self.vc[istate,0] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,1] = np.sin(phangle)/np.sqrt(2)
            self.vc[istate,2] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,3] = np.sin(phangle)/np.sqrt(2)
            self.pmod = 0
            
        elif( self.cdim == 4 and self.imod == 1 ): # 4-component split initial  (split-2 components)
            self.vc[istate,0] = np.cos(phangle)/2
            self.vc[istate,1] = np.sin(phangle)/2
            self.vc[istate,2] = np.cos(phangle)/2
            self.vc[istate,3] = np.sin(phangle)/2
            self.pmod = 1
            
        elif( self.cdim == 4 and self.imod == 2 ): # 4-component mixed initial
            self.vc[istate,0] = np.cos(phangle)
            self.vc[istate,1] = np.sin(phangle)
            self.vc[istate,2] = np.cos(phangle)
            self.vc[istate,3] = np.sin(phangle)
            self.pmod = 2
            
        elif( self.cdim == 8 and self.imod == 2): # 8-components mixed initial
            self.vc[istate,0] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,1] = np.sin(phangle)/np.sqrt(2)
            self.vc[istate,2] = -np.sin(phangle)/np.sqrt(2)
            self.vc[istate,3] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,4] = np.cos(phangle)/np.sqrt(2)
            self.vc[istate,5] = np.sin(phangle)/np.sqrt(2)
            self.vc[istate,6] = -np.sin(phangle)/np.sqrt(2)
            self.vc[istate,7] = np.cos(phangle)/np.sqrt(2)
            self.pmod = 3
        else:
            print('iniz error')
            exit()
            
    def popu(self):
        self.pop = np.zeros(self.vsize)
        self.E = 0;
        
        if ( self.pmod ==0 ):
            for i in range(self.cdim):
                self.pop += self.vc[:,i]**2
                self.E += np.dot( self.vc[:,i], np.dot(self.H, self.vc[:,i]) )
        elif (self.pmod == 1):
            self.pop = (self.vc[:,0] + self.vc[:,2])**2 + (self.vc[:,1] + self.vc[:,3])**2
            self.E = ( np.dot( self.vc[:,0]+self.vc[:,2], np.dot(self.H, self.vc[:,0]+self.vc[:,2] ) ) +
                       np.dot( self.vc[:,1]+self.vc[:,3], np.dot(self.H, self.vc[:,1]+self.vc[:,3] ) ) )
        elif (self.pmod == 2 ):
            self.pop = self.vc[:,0]*self.vc[:,2] + self.vc[:,1]*self.vc[:,3]
            self.E = ( np.dot( self.vc[:,0], np.dot(self.H, self.vc[:,2] ) ) +
                       np.dot( self.vc[:,1], np.dot(self.H, self.vc[:,3] ) ) )
        elif (self.pmod == 3 ):
            self.pop = self.vc[:,0]*self.vc[:,4] + self.vc[:,1]*self.vc[:,5] + (
                self.vc[:,2]*self.vc[:,6] + self.vc[:,3]*self.vc[:,7] )
            self.E = ( np.dot( self.vc[:,0], np.dot(self.H, self.vc[:,4] ) ) +
                       np.dot( self.vc[:,1], np.dot(self.H, self.vc[:,5] ) ) +
                       np.dot( self.vc[:,2], np.dot(self.H, self.vc[:,6] ) ) +
                       np.dot( self.vc[:,3], np.dot(self.H, self.vc[:,7] ) ) )
        else:
            print('popu loss')
            pass
        return self.pop
    
    def runmd(self):
        self.iniz()
        rcd = np.zeros((self.N, self.vsize))
        totN = np.zeros(self.N)
        totE = np.zeros(self.N)
        
        
        for i in range(self.N):
            #HVdt = np.dot(self.H, self.vc)*self.dt
            #HV1dt = np.dot(self.Hd, self.vc)*self.dt
            #HV2dt = np.dot(self.Hnd, self.vc)*self.dt
            
            if(self.rmod == 0): # 2-component equation
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                self.vc[:,1] += - np.dot(self.H, self.vc[:,0])*self.dt
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                
            elif(self.rmod == 1): # 4-component equation: Model I
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                self.vc[:,3] += - np.dot(self.H, self.vc[:,2])*self.dt/2
                
                self.vc[:,1] += - np.dot(self.H, self.vc[:,0])*self.dt                
                self.vc[:,2] += np.dot(self.H, self.vc[:,3])*self.dt
                
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                self.vc[:,3] += -np.dot(self.H, self.vc[:,2])*self.dt/2
                
            elif(self.rmod == 2): # 4-component equation: Model my version
                self.vc[:,0] -= HVdt[:,2]
                self.vc[:,1] -= HVdt[:,0]
                self.vc[:,2] += HVdt[:,3]
                self.vc[:,3] += HVdt[:,1]
                
            elif(self.rmod == 3): # 4-component equation: Model III
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/8
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/4 - np.dot(self.Hnd, self.vc[:,2])*self.dt/4
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/4
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/4 - np.dot(self.Hnd, self.vc[:,2])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/8
                
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/2 - np.dot(self.Hnd, self.vc[:,0])*self.dt/2
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/2 + np.dot(self.Hnd, self.vc[:,3])*self.dt/2
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/2
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/2 - np.dot(self.Hnd, self.vc[:,0])*self.dt/2
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/2 + np.dot(self.Hnd, self.vc[:,3])*self.dt/2
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4
                
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/8
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/4 - np.dot(self.Hnd, self.vc[:,2])*self.dt/4
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/4
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/4 - np.dot(self.Hnd, self.vc[:,2])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/8
                
            elif(self.rmod == 4): # 4-component equation: Model IV
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/8
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4 - np.dot(self.Hnd, self.vc[:,0])*self.dt/4
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/4
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4 - np.dot(self.Hnd, self.vc[:,0])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/8
                
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/2 + np.dot(self.Hnd, self.vc[:,3])*self.dt/2
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/2 - np.dot(self.Hnd, self.vc[:,2])*self.dt/2
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/2
                self.vc[:,0] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/2 + np.dot(self.Hnd, self.vc[:,3])*self.dt/2
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/2 - np.dot(self.Hnd, self.vc[:,2])*self.dt/2
                self.vc[:,1] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4
                
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/8
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4 - np.dot(self.Hnd, self.vc[:,0])*self.dt/4
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,3])*self.dt/4
                self.vc[:,2] += 0.5*np.dot(self.Hd, self.vc[:,1])*self.dt/4 + np.dot(self.Hnd, self.vc[:,1])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,0])*self.dt/4 - np.dot(self.Hnd, self.vc[:,0])*self.dt/4
                self.vc[:,3] += -0.5*np.dot(self.Hd, self.vc[:,2])*self.dt/8
                
            elif(self.rmod == 5): # 4-component equation: Model V
                print("not write")
                exit()
                
            elif(self.rmod ==6):
                print("not write")
                exit()
                
            elif(self.rmod == 7): # 8-component equation: Model X
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                self.vc[:,2] += np.dot(self.H, self.vc[:,3])*self.dt/2
                self.vc[:,5] += - np.dot(self.H, self.vc[:,4])*self.dt/2
                self.vc[:,7] += - np.dot(self.H, self.vc[:,6])*self.dt/2
                
                self.vc[:,1] += - np.dot(self.H, self.vc[:,0])*self.dt                
                self.vc[:,3] += - np.dot(self.H, self.vc[:,2])*self.dt
                self.vc[:,4] += np.dot(self.H, self.vc[:,5])*self.dt                
                self.vc[:,6] += np.dot(self.H, self.vc[:,7])*self.dt
                
                self.vc[:,0] += np.dot(self.H, self.vc[:,1])*self.dt/2
                self.vc[:,2] += np.dot(self.H, self.vc[:,3])*self.dt/2
                self.vc[:,5] += - np.dot(self.H, self.vc[:,4])*self.dt/2
                self.vc[:,7] += - np.dot(self.H, self.vc[:,6])*self.dt/2
                
            rcd[i,:] = self.popu()
            totN[i] = sum(rcd[i,:])
            totE[i] = self.E
            
            
        popdat = pd.DataFrame(rcd)
        totNdat   = pd.DataFrame(totN)
        totEdat   = pd.DataFrame(totE)
        
        popdat.to_csv('pop.dat',header=None,sep=' ')
        totNdat.to_csv('totN.dat',header=None,sep=' ')
        totEdat.to_csv('totE.dat',header=None,sep=' ')
            

idict={'dt':0.02, 'Nstep':500, 'phangle':0.5 ,'istate':0}
lamb = 0.2

Hd = np.array([[10,0,0],[0,7,0],[0,0,2]])
Hnd = lamb*np.array([[0,1,1],[1,0,1],[1,1,0]])
H = Hd+Hnd
Nsize = len(H)
Ndim = 4

mapvc = vecc(vsize=Nsize, cdim=Ndim, rmod=1, imod=2)
mapvc.getkeys(idict)
mapvc.getH(H)
mapvc.getHs(Hd,Hnd)
mapvc.runmd() 
