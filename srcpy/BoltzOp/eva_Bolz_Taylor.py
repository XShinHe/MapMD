import numpy as np
import math
H = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
O = np.eye(4)
P = np.eye(4)
beta = 100
N = 10000000
dt = beta/N
for i in range(1,N):
    P = np.dot(P,H)*beta/i
    O += (-1)**(i)*P
print(O)
