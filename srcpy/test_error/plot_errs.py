import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# note:
# csim --- complex simplectric
# cvv ---- complex velocity verlet
# qvvs --- quaternion VV & seperated term population
# qvvg --- quaternion VV & grouped term population
# qvvc --- quaternion VV & crossing term population

dt=0.01
pop_exact = pd.read_csv('pop_csim52.dat',header=None,sep='\s+')
t = pop_exact.values[:,0]*dt
pop_exact = pop_exact.values[:,1]

pop_cvv = pd.read_csv('pop_cvv22.dat',header=None,sep='\s+')
pop_cvv = np.log10(abs(pop_cvv.values[:,1]-pop_exact))
#plt.plot(t,pop_cvv,'r--',label='complex VV')

pop_qvvs = pd.read_csv('pop_qvvs22.dat',header=None,sep='\s+')
pop_qvvs = (abs( pop_qvvs.values[:,1] - pop_exact ))
plt.plot(t,pop_qvvs,'b--',label='quaternion VV seperated')

pop_qvvg = pd.read_csv('pop_qvvg22.dat',header=None,sep='\s+')
pop_qvvg = (abs( pop_qvvg.values[:,1] - pop_exact ))
plt.plot(t,pop_qvvg,'--',c=(1,0,1),label='quaternion VV grouped')

pop_qvvc = pd.read_csv('pop_qvvc22.dat',header=None,sep='\s+')
pop_qvvc = (abs( pop_qvvc.values[:,1] - pop_exact ))
plt.plot(t,pop_qvvc,'g--',label='quaternion VV crossing')

pop_ovvc = pd.read_csv('pop_ovvc22.dat',header=None,sep='\s+')
pop_ovvc = (abs( pop_ovvc.values[:,1] - pop_exact ))
plt.plot(t,pop_ovvc,'r--',label='octonion VV2 crossing')

pop_orfc = pd.read_csv('pop_orfc22.dat',header=None,sep='\s+')
pop_orfc = (abs( pop_orfc.values[:,1] - pop_exact ))
plt.plot(t,pop_orfc,'k--',label='octonion revVV crossing')

pop_csim = pd.read_csv('pop_csim22.dat',header=None,sep='\s+')
pop_csim = np.log10(abs(pop_csim.values[:,1]-pop_exact))
#plt.plot(t,pop_csim,'k--',label='complex simplectric')

plt.xlabel('t [a.u.]')
plt.ylabel(r'$Abs(Error_{population})$')
plt.legend(loc=2)
plt.savefig('errs_dt=%.3f.png'%dt)

plt.show()
