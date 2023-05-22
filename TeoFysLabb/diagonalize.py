
"""
Matrix diagonalization sample code
MW version 220325
"""

import numpy as np
from numpy import linalg as LA
from matplotlib import pyplot as plt

dim = 20 # truncated dimension of Hilbert space

def delta(m,n):
    d = 0
    if m==n: d=1
    return d

def xm(m,n):
    return (np.sqrt(n+1)*delta(m,n+1)+np.sqrt(n)*delta(m,n-1))/np.sqrt(2)

def pm(m,n):
    return 1j*(np.sqrt(n+1)*delta(m,n+1)-np.sqrt(n)*delta(m,n-1))/np.sqrt(2)

x = np.zeros(dim*dim, dtype=complex).reshape(dim, dim)
p = np.zeros(dim*dim, dtype=complex).reshape(dim, dim)
for i in range(dim):
    for j in range(dim):
        x[i,j] = xm(i,j)
        p[i,j] = pm(i,j)
x2 = x.dot(x)
p2 = p.dot(p)

# Harmonic Hamiltonian
nmax = 5
#H0 = p2/2+x2/2
H0 = p2/2+x2/2 + x2.dot(x2)
E0 = LA.eigvalsh(H0)
#print('E0=',E0[:10])
n = np.arange(dim)
E0exact = n+1/2
E02 = [0.804248812404694,2.750000000029104,5.203393244737526,8.035778658426716,11.250000000029104] #Energies from problem 2
print('E0=',E0[:nmax])

# plot eigenvalues vs n up to nmax
if nmax>dim: nmax = dim
n = np.arange(dim)
plt.figure(num=None, figsize=(12,8), dpi=80, facecolor='w', edgecolor='k')
plt.plot(n[:nmax], E0[:nmax], 'ro-', n[:nmax], E02[:nmax], 'go-')
plt.legend(['$E_n$','$E_n$ from problem 2'])
plt.xlabel('$n$')
plt.ylabel('$E_n$')
#plt.savefig('eigenvalues.pdf')
plt.show()


