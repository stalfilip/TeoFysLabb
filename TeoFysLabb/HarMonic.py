# -*- coding: utf-8 -*-
# Python simulation of a 1d harmonic oscillator
# Integrate time independent SE using the Verlet method
# Boundary conditions are found by shooting

import numpy as np
import matplotlib.pyplot as plt

hbar = 1.0 # Planck's constant (scaled)
m = 1.0    # mass (scaled)

N = 10000   # number of mesh points
a = 10.0    # range of x (-a/2 to a/2)
dx = a/N    # step length
dx2 = dx**2 # step length squared

E = 1.5    # input energy (scaled) (0.5 is the ground state energy)

# potential energy function
def V(x):
    return 0.5 * x**2 #+ x**4

# initial values and lists
x = -a/2            # initial value of position x
Even = False

if Even:
    print("Even")
    psi = 1.0
    dpsi = 0.0
else:
    print("Odd")
    psi = 0.0
    dpsi = 1.0

x_tab = []          # list to store positions for plot
psi_tab = []        # list to store wave function for plot
x_tab.append(x/a)
psi_tab.append(psi)

c = 2.0*m/hbar**2 # constant in Schrödinger equation

for i in range(N) :
    d2psi = c*(V(x)-E)*psi
    d2psinew = c*(V(x+dx)-E)*psi
    psi += dpsi*dx + 0.5*d2psi*dx2
    dpsi += 0.5*(d2psi+d2psinew)*dx
    x += dx
    x_tab.append(x/a)
    psi_tab.append(psi)

print('E = ', E, ', psi(x=a/2) = ', psi)

plt.close()
plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
plt.plot(x_tab, psi_tab, linewidth=1, color='red')
plt.xlabel('x')
plt.ylabel('$\psi$')
plt.show()