# -*- coding: utf-8 -*-
# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Boundary conditions are found by shooting
# MW 220325

import numpy as np
import matplotlib.pyplot as plt

a=1.0e-9           # well width a=1 nm
hbar=1.0545718e-34 # Plancks constant
m=9.10938356e-31   # electron mass
e=1.60217662e-19   # electron charge=-e
c=2.0*m/hbar**2    # constant in Schrödinger equation
N=10               # number of mesh points
dx=a/N             # step length
dx2=dx**2          # step length squared

EeV=(hbar*np.pi/a)**2/(2.0*m)/e
#EeV = 0.0          # input energy in eV: test 0.3 , 0.4 , 0.3760 , 1.5
E = EeV*e          # input energy in J

#print('Exact solution for infinite box potential')
print('E1=',(hbar*np.pi/a)**2/(2*m)/e,'eV')
#print('E2=',(hbar*2*np.pi/a)**2/(2*m)/e,'eV')

# potential energy function
def V(x):
    y = 0.0
    #y = x**2/2 # harmonic oscillator
    #y = x**2/2 + x**4 # anharmonic oscillator
    return y

# initial values and lists
x = 0               # initial value of position x
psi = 0.0           # wave function at initial position
dpsi = 1.0          # derivative of wave function at initial position
x_tab = []          # list to store positions for plot
psi_tab = []        # list to store wave function for plot
x_tab.append(x/a)
psi_tab.append(psi)

for i in range(N) :
    d2psi = c*(V(x)-E)*psi
    d2psinew = c*(V(x+dx)-E)*psi
    psi += dpsi*dx + 0.5*d2psi*dx2
    dpsi += 0.5*(d2psi+d2psinew)*dx
    x += dx
    x_tab.append(x/a)
    psi_tab.append(psi)

print('E=',EeV,'eV , psi(x=a)=',psi)

plt.close()
plt.figure(num=None, figsize=(8,8), dpi=80, facecolor='w', edgecolor='k')
plt.plot(x_tab, psi_tab, linewidth=1, color='red')
#plt.xlim(0, 1)
#limit=1.e-9
#plt.ylim(0, limit)
#plt.ylim(-limit, limit)
#plt.autoscale(False)
plt.xlabel('x/a')
plt.ylabel('$\psi$')
#plt.savefig('psi.pdf')
plt.show()

