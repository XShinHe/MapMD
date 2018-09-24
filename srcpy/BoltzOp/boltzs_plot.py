import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a=pd.read_csv("boltzs.txt")
rcd=a.values
bx=1
ex=180
plt.plot(rcd[0,bx:ex],rcd[1,bx:ex],'b--',marker='o',markersize=3,label='Dynamics I')
plt.plot(rcd[0,bx:ex],rcd[2,bx:ex],'g--',marker='s',markersize=3,label='Dynamics II')
plt.plot(rcd[0,bx:ex],rcd[3,bx:ex],'r--',marker='x',markersize=3,label='Taylor series')
plt.title("Comparision of evaluation Boltzmann Operator \n(calculation loop in N=10000 step)")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$Tr[e^{-\beta H}]$")
plt.legend(loc=1)
plt.savefig('boltzs.png')
plt.show()
