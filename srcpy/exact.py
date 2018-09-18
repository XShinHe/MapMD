import numpy as np
M = np.array([[0.06,0.22,-0.22,0.08],[0.22,1.12,-0.08,0.25],[-0.22,-0.08,1.12,-0.25],[0.08,0.25,-0.25,2.04]])
print('Hamiltonian \n',M)
eige,eigv=np.linalg.eigh(M)
print('eigenvalues \n',eige)
print('eigenvectors \n',eigv)
