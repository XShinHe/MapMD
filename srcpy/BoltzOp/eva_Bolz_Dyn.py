import numpy as np
H = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
O = np.eye(4)
beta = 0.01
N = 100000
dt = beta/N
for i in range(N):
    O += -np.dot(H,O)*dt
print(O)
