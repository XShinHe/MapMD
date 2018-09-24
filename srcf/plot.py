import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

a = pd.read_csv('coo.dat',sep='\s+',header=None)
a = a.values.T
plt.plot(a[1],a[2],'r.')
plt.savefig('phasespace1.png')
plt.show()

